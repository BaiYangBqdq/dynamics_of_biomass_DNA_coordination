% this script describe the hns dnaA initial model + adder division model
% the core of this model is OriC/determine adder at each division

close all ; clear all ;
tic;
%% global parameters
TT = 1000 ; % total calculation time (min)
dt = 0.1 ; % calculation time step (min)
dtR = 1 ; % time step for recording (min)
RTlist = dtR : dtR : TT ; % time list (min)
% Ctime = 250; % changing time of hns expression (min)


CAcc = 10 ; % critical cytoplasmic dnaA concentration to initiate DNA replication 
alpha0 = 1;
datAn =100 ; % dnaA binding site on datA
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

CNi = 1000; % initial cell number
CNmax = 10000; % maxi cell number
CNdilut = 0.8*CNmax ; % dilute the cells at this cell number
dilutRatio = 10 ; 
Tdilut = nan(1,20) ; % diluting ratio % dilution time

Tlist = dt:dt:TT;
% alpha0_list = nan(size(Tlist));
% alpha0_list(1:Ctime/dt) = alpha0;
% alpha0_list(Ctime/dt:(Ctime+2*DT)/dt) = alpha0:-(0.8*alpha0/((Ctime+2*DT)/dt-Ctime/dt)):0.2*alpha0;
% alpha0_list((Ctime+2*DT)/dt:end) = 0.2*alpha0;
alpha0_list = alpha0*ones(size(Tlist));

% figure;
% plot(Tlist,alpha0_list);
% xlabel('time');
% ylabel('\alpha_0')
% drawnow;


%% recording vectors

Rec.A_conc_mean = nan(1,length(RTlist)) ; % mean cytoplasmic dnaA protein
Rec.A_conc_std = nan(1,length(RTlist)) ; % std cytoplasmic dnaA protein
Rec.V_mean = nan(1,length(RTlist)) ; % mean cell volume ;
Rec.V_std = nan(1,length(RTlist)) ; % std cell volume ;
Rec.DNA = nan(1,length(RTlist)) ; % cellular DNA mean
Rec.OriC = nan(1,length(RTlist)) ; % total oric
Rec.CN = nan(1,length(RTlist)) ; % total cell number
Rec.CNs = nan(1,length(RTlist)) ; % total cell number
Rec.V = nan(1,length(RTlist)) ; % total mass
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

% Cell = repmat(C,CNmax, 1) ;
%% initiation 

Cell(1).DNAcopy = 1*ones(1,DNAlength);
Cell(1).DNA = sum(Cell(1).DNAcopy)/DNAlength ;
Cell(1).OriC = Cell(1).DNAcopy(1);
Cell(1).Nd = Cell(1).DNAcopy(dTdatA/dt+1); 

Cell(1).nrep = 1 ; % number of DNA replication. starting from 1, not 0
Cell(1).ndiv = 1; % number of cell division
Cell(1).RpO = 1 ; % replication initiation operator, 1 for open, 0 for closed
Cell(1).V = 1 ; % initial cell mass
Cell(1).A = 1 ; % initial dnaA protein
Cell(1).H = alphaH0 * Cell(1).V;
Cell(1).Abox = max((Cell(1).Nd*datAn + Cell(1).DNA*Chns - Cell(1).H),0);
Cell(1).Af = max((Cell(1).A - Cell(1).Abox),0) ;
Cell(1).lambda = lambda0;
Cell(1).alpha = alpha0; 

Cell(1).Prep = NaN(1,ceil(TT/DT)+10)  ; % position of DNA replication  
Cell(1).RpO = 1 ; 
Cell(1).ndiv = 1;

Cell(1).rpt = NaN(1,ceil(TT/DT)+10) ;
Cell(1).divt =  NaN(1,ceil(TT/DT)+10)  ; % division time
Cell(1).Prep = NaN(1,ceil(TT/DT)+10)  ; % position of DNA replication  
Cell(1).Tdiv = [inf];


%% main loop
CN = 1 ;
ndilut = 0 ; % number of dilution times

for t = 1 : length(Tlist)
    
    alpha0 = alpha0_list(t);
    
    for i = 1 : CN
    Cell(i).DNA = sum(Cell(i).DNAcopy)/DNAlength ;
    Cell(i).Nd = Cell(i).DNAcopy(dTdatA/dt+1); % copy number of datA
    Cell(i).A = Cell(i).A + Cell(i).alpha*Cell(i).V*dt;
    Cell(i).V = Cell(i).V + Cell(i).lambda*Cell(i).V*dt ;
    Cell(i).H = alphaH0*Cell(i).V;
    Cell(i).Abox = max((Cell(i).Nd*datAn + Cell(i).DNA*Chns - Cell(i).H),0);
    Cell(i).Af = max(Cell(i).A-Cell(i).Abox,0);
    
    
    % check new replication initiation condition
    if Cell(i).Af/Cell(i).V > CAcc && Cell(i).RpO == 1
        Cell(i).nrep = Cell(i).nrep +1; % replication number +1
        Cell(i).Prep(Cell(i).nrep) = 1 ; % generate new replication position
        Cell(i).rpt(Cell(i).nrep) = t*dt; % replication time +dt
        Cell(i).RpO = 0; % stop new replicaiton initiation
        Cell(i).Tdiv(Cell(i).nrep) = Cell(i).rpt(Cell(i).nrep) + Tcycle;
        
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
    if any(Tlist(t) >= Cell(i).Tdiv) 
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
    
    if CN>CNdilut
       sampleindex = randsample(CN,round(CN/dilutRatio));
       Cellcopy = Cell ;
       Cell = repmat(C,CNmax, 1) ;
       Cell(1:round(CN/dilutRatio)) = Cellcopy(sampleindex);
       CN = round(CN/dilutRatio) ; % 
       ndilut = ndilut +1 ; % number of dilution times
       disp(['dilution time = ' num2str(t*dt) '/' num2str(TT) ', CN = ' num2str(CN)]); % display  time
       Tdilut(ndilut) = t*dt ; % dilution time
       clear Cellcopy
   end
    % recording
    if rem(t*dt,dtR)==0 
        Rec.A_conc_mean(t*dt/dtR) = mean([Cell(1:CN).A]./[Cell(1:CN).V]) ; % mean cytoplasmic dnaA protein
        Rec.A_conc_std(t*dt/dtR) = std([Cell(1:CN).A]./[Cell(1:CN).V]) ; % std cytoplasmic dnaA protein
        Rec.V_mean(t*dt/dtR) = mean([Cell(1:CN).V])  ; % mean cell volume ;
        Rec.V_std(t*dt/dtR) = std([Cell(1:CN).V])  ; % std cell volume ;

        Rec.OriC(t*dt/dtR) = nansum([Cell(1:CN).OriC]) ; 
        Rec.DNA(t*dt/dtR) = nansum([Cell(1:CN).DNA]) ; 
        Rec.V(t*dt/dtR) = nansum([Cell(1:CN).V]) ; 
        
        Rec.CN(t*dt/dtR) = CN*dilutRatio^ndilut ; % total cell number
        Rec.CNs(t*dt/dtR) = CN ; % cell number in the testing sample
    end
    
    
end

toc;
%% plot time course data

hp=figure(); % plot cellualr proteins
semilogy(RTlist,exp(lambda0*RTlist),'linewidth',2); hold all
semilogy(RTlist,Rec.OriC,'linewidth',2); hold all
semilogy(RTlist,Rec.DNA,'linewidth',2); hold all
semilogy(RTlist,Rec.CN,'linewidth',2); hold all

xlabel('time');legend('mass','oriC','DNA','cell number','location','NorthWest'); hold off
 drawnow

 
 %% plot iDR
hiDR=figure();
yaxis = [0 4] ;
winds = 30; % smooth winds, smooth over these many points
SCN = smoothts(Rec.CN, 'b', winds) ; % smooth CN
SiCDR = 60*(log(SCN(2:end))-log(SCN(1:end-1)))/dtR; % iCDR
SOriC = smoothts(Rec.OriC, 'b', winds) ; % smooth mean OriC
SiODR = 60*(log(SOriC(2:end))-log(SOriC(1:end-1)))/dtR; % iODR
SV = smoothts(Rec.V, 'b', winds) ; % smooth mean V
SiVDR = 60*(log(SV(2:end))-log(SV(1:end-1)))/dtR; % iVDR
SG = smoothts(Rec.DNA, 'b', winds) ; % smooth mean DNA
SiGDR = 60*(log(SG(2:end))-log(SG(1:end-1)))/dtR; % iGDR

subplot(4,1,1)
plot(RTlist(1:end-1),SiVDR,'linewidth',2); hold all; title('iVDR')
subplot(4,1,2)
plot(RTlist(1:end-1),SiCDR,'linewidth',2); hold all; title('iCDR')
subplot(4,1,3)
plot(RTlist(1:end-1),SiODR,'linewidth',2); hold all; title('iODR')
subplot(4,1,4)
plot(RTlist(1:end-1),SiGDR,'linewidth',2); hold all; title('iGDR')
xlabel('Time (min)'); 
 

%%
CS = Cell(1:CN);
save('Cell_steady_state.mat','CS');

