% this script describe the hns dnaA initial model + adder division model
% the core of this model is OriC/determine adder at each division

close all ; clear all ;
tic;
%% global parameters
TT = 1000 ; % total calculation time (min)
dt = 0.1 ; % calculation time step (min)
dtR = 1 ; % time step for recording (min)
RTlist = dtR : dtR : TT ; % time list (min)
Ctime = 500; % changing time of hns expression (min)


CAcc = 10 ; % critical cytoplasmic dnaA concentration to initiate DNA replication 
alpha0 = 5;
datAn = 100 ; % dnaA binding site on datA
Chns = 300 ; % number of hns binding to each chromosome
alphaH0 = 200;
Tseq = 10;
dTdatA = 10 ; % time delay between OriC and datA (min)

DT = 23; 
Cperiod = 42;
DNAlength = ceil(Cperiod/dt) ; % DNA length is set to fix V fork = 1 ;
lambda0 = log(2)/DT ; % cell mass growth rate (min^-1)


Tlist = dt:dt:TT;
% alpha0_list = nan(size(Tlist));
% alpha0_list(1:Ctime/dt) = alpha0;
% alpha0_list(Ctime/dt:end) = 0;


% figure;
% plot(Tlist/60,alpha0_list);
% xlabel('time (hr)');
% ylabel('\alpha_0')
% drawnow;

%% initiation 

DNAcopy = 1*ones(1,DNAlength);
DNA = sum(DNAcopy)/DNAlength ;
OriC = DNAcopy(1);
Nd = DNAcopy(dTdatA/dt+1); 

nrep = 1 ; % number of DNA replication. starting from 1, not 0
ndiv = 1; % number of cell division

V = 1 ; % initial cell mass
A = 1 ; % initial dnaA protein
H = alphaH0 * V;
Abox = max((Nd*datAn + DNA*Chns - H),0);
Af = max((A - Abox),0) ;
lambda = lambda0;
alpha = alpha0; 

rpt = NaN(1,ceil(TT/DT)+10) ;
Prep = 1;
RpO = 1;

%% main loop
for t = 1 : length(Tlist)
    
%     alpha = alpha0_list(t);
    alpha = alpha0;
    
    OriC = DNAcopy(1);
    DNA = sum(DNAcopy)/DNAlength ;
    Nd = DNAcopy(dTdatA/dt+1); % copy number of datA
    A = A + alpha*V*dt;
    V = V + lambda*V*dt ;
    H = alphaH0*V;
    Abox = max((Nd*datAn + DNA*Chns - H),0);
    Af = max(A-Abox,0);
    
    % check new replication initiation condition
    if Af/V > CAcc && RpO == 1
%     if Af/V > CAcc
        nrep = nrep +1; % replication number +1
        Prep(nrep) = 1 ; % generate new replication position
        rpt(nrep) = t*dt; % replication time +dt
        RpO = 0; % stop new replicaiton initiation
    end
    
    if t == ceil((rpt(nrep) + Tseq)/dt); RpO = 1; end
    
    
    % DNA replication
    for f = 1 : length(Prep) 
        if Prep(f) <= Cperiod/dt
           DNAcopy(Prep(f)) = DNAcopy(Prep(f))*2 ; % all DNA replicates
            Prep(f) = Prep(f) +1 ; % replication position moves
        else if Prep(f) > Cperiod/dt % DNA replication terminates
            Prep(f) = NaN ;   
            end
        end
    end
    

    % recording
    if rem(t*dt,dtR)==0 
        Rec_A(t*dt/dtR) = A ; % 
        Rec_Af(t*dt/dtR) = Af ; % 
        Rec_V(t*dt/dtR) = V  ; 
        Rec_OriC(t*dt/dtR) = OriC ; 
        Rec_DNA(t*dt/dtR) = DNA ; 
    end
    
    
    list_A(t) = A ; % 
    list_Af(t) = Af ; % 
    list_V(t) = V ; % 
    list_OriC(t) = OriC ; % 
    
end

toc;
%% plot time course data

ht=figure(); 
semilogy(RTlist/60,Rec_V/(Rec_V(Ctime/dtR)),'linewidth',2); hold all
semilogy(RTlist/60,Rec_OriC/(Rec_OriC(Ctime/dtR)),'linewidth',2); hold all
semilogy([Ctime Ctime]/60,[10^-10 10^5],'k--','linewidth',2); hold all
xlim([0 15]); ylim([10^-10 10^5])
xlabel('Time (hr)'); ylabel('Rel. value');
legend('mass','oriC','location','NorthWest'); hold off
set(gca,'fontsize',14);
drawnow


%%
h=figure('position',[100 100 500 600]); 
% subplot(3,1,1)
% plot((Tlist(3000:7000)-Ctime)/60,alpha0_list(3000:7000),'linewidth',2); hold all
% ylabel('\alpha_A'); xlim([-3 3]);
% set(gca,'fontsize',14);
subplot(2,1,1)
plot((Tlist(3000:7000)-Ctime)/60,list_Af(3000:7000)./list_V(3000:7000),'linewidth',2); hold all
plot((Tlist(3000:7000)-Ctime)/60, 10*ones(size((Tlist(3000:7000)-Ctime)/60)),'k--','linewidth',2);
ylabel('free dnaA / mass');xlim([-3 3]);
set(gca,'fontsize',14);
subplot(2,1,2)
plot((Tlist(3000:7000)-Ctime)/60,list_OriC(3000:7000)./list_V(3000:7000),'linewidth',2); hold all
ylabel('oriC / mass');xlim([-3 3]);

xlabel('Time (hr)');
set(gca,'fontsize',14);
drawnow
