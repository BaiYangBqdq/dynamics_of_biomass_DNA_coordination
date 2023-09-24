% this script describe the hns dnaA initial model + adder division model
% the core of this model is OriC/determine adder at each division

close all ; clear all ;
tic;
%% global parameters
TT = 240 ; % total calculation time (min)
dt = 0.1 ; % calculation time step (min)
dtR = 1 ; % time step for recording (min)
RTlist = dtR : dtR : TT ; % time list (min)
Ctime = 60; % changing time of hns expression (min)


CAcc = 10 ; % critical cytoplasmic dnaA concentration to initiate DNA replication 
alpha0 = 1;
datAn = 100 ; % dnaA binding site on datA
Chns = 300 ; % number of hns binding to each chromosome
alphaH0 = 50;
Tseq = 10;
dTdatA = 0 ; % time delay between OriC and datA (min)

DT = 23; 
Cperiod = 42;
Tcycle = 63;
DNAlength = ceil(Cperiod/dt) ; % DNA length is set to fix V fork = 1 ;
lambda0 = log(2)/DT ; % cell mass growth rate (min^-1)

epsilon_A = 0.1;
epsilon_V = 0.1;

CNmax = 200000; % maxi cell number

Tlist = dt:dt:TT;
alpha0_list = nan(size(Tlist));
alpha0_list(1:Ctime/dt) = alpha0;
alpha0_list(Ctime/dt:(Ctime+2*DT)/dt) = alpha0:-(0.8*alpha0/((Ctime+2*DT)/dt-Ctime/dt)):0.2*alpha0;
alpha0_list((Ctime+2*DT)/dt:end) = 0.2*alpha0;


figure;
plot((Tlist-Ctime)/60,alpha0_list,'linewidth',2); hold all
plot([0 0],[0 1.5],'k--','linewidth',2); hold all
xlabel('time (hr)');
ylabel('\alpha_0');
ylim([0 15]);
drawnow;


%% recording vectors

Rec.A_conc_mean = nan(1,length(RTlist)) ; % mean cytoplasmic dnaA protein
Rec.A_conc_std = nan(1,length(RTlist)) ; % std cytoplasmic dnaA protein
Rec.DNA = nan(1,length(RTlist)) ; % cellular DNA mean
Rec.OriC = nan(1,length(RTlist)) ; % total oric
Rec.V = nan(1,length(RTlist)) ; % total mass
Rec.A = nan(1,length(RTlist)) ; % total A protein
Rec.Abox = nan(1,length(RTlist)) ; % total A box
Rec.CN = nan(1,length(RTlist)) ; % total cell number


%% single cell phisiological properties structure

C.A = nan ; % cellular dnaA protein ; 
C.Af = nan ; % cytoplasmic dnaA protein
C.H = nan ; % cytoplasmic hns protein
C.Abox = nan;
C.V = nan ; % cell volume ;

C.DNAcopy = ones(1,DNAlength); % DNA copy;
C.OriC = C.DNAcopy(1); % cellular OriC
C.Nd = C.DNAcopy(dTdatA/dt+1); % copy number of datA
C.DNA = sum(C.DNAcopy)/DNAlength ; % cellular DNA amount

C.nrep = nan ; % number of DNA replication. starting from 1, not 0
C.ndiv = nan;
C.RpO = nan ; % replication initiation operator, 1 for open, 0 for closed
C.lambda = nan; 
C.alpha = nan; 
C.RpO = nan ; % replication initiation operator, 1 for open, 0 for closed

C.rpt = NaN(1,ceil(TT/DT)+10) ;
C.divt =  NaN(1,ceil(TT/DT)+10)  ; % division time
C.Prep = NaN(1,ceil(TT/DT)+10)  ; % position of DNA replication  
C.Tdiv = [inf];

Cell = repmat(C,CNmax, 1) ;
%% initiation 

load Cell_steady_state.mat;

CN=length(CS);
Cell(1:CN)=CS;
for i = 1:CN
    Cell(i).rpt = Cell(i).rpt - 1000;
    Cell(i).divt = Cell(i).divt - 1000;
    Cell(i).Tdiv = Cell(i).Tdiv - 1000;
end



%% main loop
for t = 1 : length(Tlist)
    
    alpha0 = alpha0_list(t);
    
    if t==Ctime/dt
        disp('dnaA shutdown');
        for i=1:CN
        Cell(i).alpha = alpha0*(epsilon_A * randn(1)+1);
        end
    end
    
    for i = 1 : CN
    Cell(i).DNA = sum(Cell(i).DNAcopy)/DNAlength ;
    Cell(i).OriC = Cell(i).DNAcopy(1); % copy number of OriC
    Cell(i).Nd = Cell(i).DNAcopy(dTdatA/dt+1); % copy number of datA
    Cell(i).A = Cell(i).A + Cell(i).alpha*Cell(i).V*dt;
    Cell(i).V = Cell(i).V + Cell(i).lambda*Cell(i).V*dt ;
    Cell(i).H = alphaH0*Cell(i).V;
    Cell(i).Abox = max((Cell(i).Nd*datAn + Cell(i).DNA*Chns - Cell(i).H),0);
    Cell(i).Af = max(Cell(i).A-Cell(i).Abox,0);
    
    
        
    % check new replication initiation condition
    if Cell(i).Af/Cell(i).V > CAcc  && Cell(i).RpO == 1
%     if Cell(i).Af/Cell(i).V > CAcc 
        Cell(i).nrep = Cell(i).nrep +1; % replication number +1
        Cell(i).Prep(Cell(i).nrep) = 1 ; % generate new replication position
        Cell(i).rpt(Cell(i).nrep) = t*dt; % replication time +dt
        Cell(i).RpO = 0; % stop new replicaiton initiation
        Cell(i).Tdiv(Cell(i).nrep) = Cell(i).rpt(Cell(i).nrep) + Tcycle;
%         disp('DNA rep init')
    end
    
    
    % DNA replication
    for f = 1 : length(Cell(i).Prep) 
        if Cell(i).Prep(f) <= Cperiod/dt
           Cell(i).DNAcopy(Cell(i).Prep(f)) = Cell(i).DNAcopy(Cell(i).Prep(f))*2 ; % all DNA replicates
            Cell(i).Prep(f) = Cell(i).Prep(f) +1 ; % replication position moves
        else if Cell(i).Prep(f) > Cperiod/dt % DNA replication terminates
            Cell(i).Prep(f) = NaN ;   
            end
        end
    end

    
    if t == ceil((Cell(i).rpt(Cell(i).nrep) + Tseq)/dt); Cell(i).RpO = 1; end
    
    % Division
    if any(Tlist(t) >= Cell(i).Tdiv) && Cell(i).DNAcopy(end)>1
        Cell(i).ndiv = Cell(i).ndiv +1;
        Cell(i).divt(Cell(i).ndiv) = t*dt;
        Cell(i).Tdiv(Tlist(t)>=Cell(i).Tdiv) = inf;
        
        Cell(i).V = Cell(i).V/2;
        Cell(i).A = Cell(i).A/2;
        Cell(i).Abox = Cell(i).Abox/2;
        Cell(i).Af = Cell(i).Af/2;
        Cell(i).DNAcopy = Cell(i).DNAcopy/2;
        Cell(i).H = Cell(i).H/2;
        
        Cell(i).OriC = Cell(i).DNAcopy(1);
        Cell(i).DNA = sum(Cell(i).DNAcopy)/DNAlength ;
        Cell(i).Nd = Cell(CN).DNAcopy(dTdatA/dt+1); 
        
        Cell(i).lambda = lambda0*(epsilon_V * randn(1)+1);
        Cell(i).alpha = alpha0*(epsilon_A * randn(1)+1);
        
        Cell(i).Prep = Cell(i).Prep  ;
        Cell(i).RpO = Cell(i).RpO ; 
        Cell(i).nrep = Cell(i).nrep ; 
        Cell(i).rpt = Cell(i).rpt;
        
        
        % add new cell
        CN = CN+1 ; % cell number +1
        Cell(CN) = Cell(i);
        Cell(CN).lambda = lambda0*(epsilon_V * randn(1)+1);
        Cell(CN).alpha = alpha0*(epsilon_A * randn(1)+1);
    end
    end
    

    % recording
    if rem(t*dt,dtR)==0 
        Rec.A_conc_mean(t*dt/dtR) = mean([Cell(1:CN).A]./[Cell(1:CN).V]) ; % mean cytoplasmic dnaA protein
        Rec.A_conc_std(t*dt/dtR) = std([Cell(1:CN).A]./[Cell(1:CN).V]) ; % std cytoplasmic dnaA protein
        Rec.A(t*dt/dtR) = nansum([Cell(1:CN).A]) ; 
        Rec.Abox(t*dt/dtR) = nansum([Cell(1:CN).Abox]) ; 
        Rec.V(t*dt/dtR) = nansum([Cell(1:CN).V]) ; 
        Rec.OriC(t*dt/dtR) = nansum([Cell(1:CN).OriC]) ; 
        Rec.DNA(t*dt/dtR) = nansum([Cell(1:CN).DNA]) ; 
        
        Rec.CN(t*dt/dtR) = CN; % total cell number
        disp(['rec time = ' num2str(t*dt) '/' num2str(TT) ', CN = ' num2str(CN) ', OriC = ' num2str(Rec.OriC(t*dt/dtR)) ', A = ' num2str(Rec.A(t*dt/dtR)) ', Abox = ' num2str(Rec.Abox(t*dt/dtR))]); % display  time
    end
    
    
end

toc;
%% plot time course data

hp=figure(); % plot cellualr proteins
yyaxis left
semilogy(RTlist/60,Rec.V,'linewidth',2); hold all
semilogy(RTlist/60,Rec.OriC*4,'linewidth',2); hold all
legend('mass','oriC','DNA','cell number','location','NorthWest'); hold off

yyaxis right
plot(Tlist/60,alpha0_list,'linewidth',2); 
xlabel('Time (hr)');

drawnow


%% plot iDR
hiDR=figure();
yaxis = [0 4] ;
winds = 2; % smooth winds, smooth over these many points
SCN = smoothts(Rec.CN, 'b', winds) ; % smooth CN
SiCDR = 60*(log(SCN(2:end))-log(SCN(1:end-1)))/dtR; % iCDR
SOriC = smoothts(Rec.OriC, 'b', winds) ; % smooth mean OriC
SiODR = 60*(log(SOriC(2:end))-log(SOriC(1:end-1)))/dtR; % iODR
SV = smoothts(Rec.V, 'b', winds) ; % smooth mean V
SiVDR = 60*(log(SV(2:end))-log(SV(1:end-1)))/dtR; % iVDR
SG = smoothts(Rec.DNA, 'b', winds) ; % smooth mean DNA
SiGDR = 60*(log(SG(2:end))-log(SG(1:end-1)))/dtR; % iGDR


plot((RTlist(2:end)-Ctime)/60,SiVDR,'linewidth',2);hold all
plot((RTlist(2:end)-Ctime)/60,SiODR,'linewidth',2);
plot([Ctime Ctime]-Ctime,yaxis,'k--'); hold off;
ylabel('iDR'); ylim(yaxis);
legend('iVDR','iODR');
xlabel('Time (min)'); 
drawnow
%%
save('Population_titration_extrusion_dnaA_shutdown.mat','-v7.3');

