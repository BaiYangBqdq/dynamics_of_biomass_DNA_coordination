% titration switch

% close all ;
clear all ;

tic;
dnaA_list = 2.^[0:0.2:6];
mean_IM = nan(size(dnaA_list));

%% global parameters
for id = 1:length(dnaA_list)

DT = 0.5;     %(h)
TT = DT*300 ; % total calculation time (h)
dt = 0.01 ;    % calculation time step (h)
dtR = 0.1;
Cperiod = 0.7; 

dTdatA = 0 ; % time delay between OriC and datA (h)
dTdars1 = 0.1 ; % time delay between OriC and datA (h)
dTdars2 = 0.2 ; % time delay between OriC and datA (h)
dTdnaA = 0 ; % time delay between OriC and dnaA (min)

% repR = 1/Cperiod ; % DNA replication rate = 1 / C period
DNAlength = Cperiod/dt ; % DNA length is set to fix V fork = 1 ;
lambda0 = log(2)/DT ; % cell mass growth rate (h^-1)

D_total = 400;
KPD = 400;
KD = 50;

alphal = 7500;
alphad1 = 1000;
alphad2 = 6430;
beta_datA = 6000;
beta_rida = 50000;
alphaH0 = 500;

% alphal = 0;
% alphad1 = 0;
% alphad2 = 0;
% beta_datA = 0;
% beta_rida = 0;


lambda = lambda0;
tseq = 15/60; % hr

phiP0 = dnaA_list(id)/50 * 10^-3 * 10^6 ;
CAcc = 20  ; % critical cytoplasmic dnaA concentration to initiate DNA replication 
n = 5;
% KPD = 200/D_total;
Chns = 0;

%% initiation 
Tlist = (1 : TT/dt)*dt ; % unit [min]
Vlist = nan(1,TT/dt) ;
CDATPlist = nan(1,TT/dt) ;
DNAcopy = ones(TT/dt,DNAlength); % starting DNAcopy = 1 ;
DNAcopy(1,1) = 1 ;
DNAlist = NaN(TT/dt,1); % total DNA
OriClist = NaN(TT/dt,1); % total DNA
rptlist = NaN(TT/dt,ceil(TT/DT)+10) ; % DNA replication time
tertlist = NaN(TT/dt,ceil(TT/DT)+10) ; % DNA termination time
Preplist = NaN(TT/dt,ceil(TT/DT)+10) ; % position of DNA replication  
NAlist = NaN(TT/dt,1); % copy number of dnaA
Ndlist = NaN(TT/dt,1); % copy number of datA
nreplist = NaN(TT/dt,1); % number of DNA replication. starting from 1, not 0
RpOlist = NaN(TT/dt,1); % replication initiation operator, 1 for open, 0 for closed
Adlist = NaN(TT/dt,1); % the total capacity that datA could absorb dnaA


NA = 1 ; % copy number of dnaA
Ndat = 1 ; % copy number of datA
nrep = 1 ; % number of DNA replication. starting from 1, not 0
RpO = 1 ; % replication initiation operator, 1 for open, 0 for closed
V = 1 ; % initial cell mass
NDATP = 100;
% fATP = 0.8;
NT = D_total*V;
NDATPf = 0;
DNA = sum(DNAcopy(1,:))/DNAlength ;
H = alphaH0 * V;
Abox = max((DNA*Chns - H),0);

rpt = NaN(1,ceil(TT/DT)+10) ; % DNA replication time
tert = NaN(1,ceil(TT/DT)+10) ; % DNA termination time
Prep = NaN(1,ceil(TT/DT)+10) ; % position of DNA replication  
Vi =  NaN(1,ceil(TT/DT)+10) ; % initiation mass
Vt =  NaN(1,ceil(TT/DT)+10) ; % termination mass
S = 1;
for t = 1 : Tlist(end)/dt
    
    V = V + lambda*V*dt ;
    H = alphaH0 * V;
    
    DNA = sum(DNAcopy(t,:))/DNAlength ;
    No = DNAcopy(t,dTdnaA/dt+1); % copy number of dnaA
    Ndat = DNAcopy(t,dTdatA/dt+1); % copy number of datA
    Nd1 = DNAcopy(t,dTdars1/dt+1); % copy number of dars1
    Nd2 = DNAcopy(t,dTdars2/dt+1); % copy number of dars2
      
    CDATP = NDATP/V;
    CDADP = (NT-NDATP)/V;
    
    Abox = max((DNA*Chns - H),0);
    NDATPf = max(NDATP - Abox,0) ;
    NTf = max(NT - Abox,0) ;
    CDf = NTf/V;
    CDATPf = NDATPf/V;   

    NT = NT + dt*(phiP0*lambda*V/(1+(CDf/KPD)^n));
    NDATP = NDATP + dt * ((alphal*V + alphad1*Nd1 + alphad2*Nd2) * CDADP/(KD+CDADP) - (beta_datA + beta_rida)*No * CDATPf/(KD+CDATPf) + phiP0*lambda*V/(1+(CDf/KPD)^n) );

    DNAcopy(t+1,:) = DNAcopy(t,:) ;
    % check new replication initiation condition
%     if CATPf > CAcc && RpO == 1
    if NDATP/NT > 0.8 && RpO == 1
        nrep = nrep +1; % replication number +1
        Prep(nrep) = 1 ; % generate new replication position
        rpt(nrep) = t*dt; % replication time +dt
        RpO = 0; % stop new replicaiton initiation
        Vi(nrep) = V/DNAcopy(t,1); % note initiation mass
%         lambda = lambda0 * (epsilon * randn(1,1)+1);
    end
    % DNA replication
    for f = 1 : length(Prep) 
        if Prep(f) <= Cperiod/dt
            DNAcopy(t+1,Prep(f)) = DNAcopy(t+1,Prep(f))*2 ; % all DNA replicates
            Prep(f) = Prep(f) +1 ; % replication position moves
        else if Prep(f) > Cperiod/dt % DNA replication terminates
            Prep(f) = NaN ;   
            tert(nrep) = t*dt ; % DNA termination time
            Vt(nrep) = V/DNAcopy(t,end-1) ; % termination mass
            end
        end
    end
    % set replication initiation operator
    if t == ceil((rpt(nrep) + tseq)/dt); RpO = 1; end
    
    
    % save data in vectors

    OriC = DNAcopy(t,1) ;
%     Af = NDATPf;
    
%     Vlist(t) = V ;
%     OriClist(t) = DNAcopy(t,1) ;
%     DNAlist(t) = DNA;
%     NTlist(t) = NT;
%     fATPlist(t) = fATP;
%     nreplist(t) = nrep;
%     RpOlist(t) = RpO; 
%     Df_list(t) = CDf;
%     DATPf_list(t) = CDATPf;
    
    
%     % recording
%     if rem(t*dt,dtR)==0 
%         Rec_fATP(round(t*dt/dtR)) = fATP ; % 
%         Rec_Af(round(t*dt/dtR))  = Af ; % 
%         Rec_V(round(t*dt/dtR))  = V  ; 
%         Rec_OriC(round(t*dt/dtR))  = OriC ; 
%         Rec_DNA(round(t*dt/dtR))  = DNA ; 
%     end
    
    CDATPflist(t) = CDATPf;
    list_V(t) = V ; % 
    list_OriC(t) = OriC ; % 
    list_t(t) = t*dt;
    
end

MpO = list_V./list_OriC;
% figure;
% plot(list_t,MpO);hold all
% % plot(list_t,CDATPflist);
% drawnow;

try
    rpt(isnan(rpt))=[];
    % G_conc = Rec_DNA./Rec_V;
    
    % mean_G_conc(id) = mean(G_conc(rpt(end-1)/dtR:rpt(end)/dtR));
    mean_IM(id) = mean(MpO(round(rpt(end-10:end)/dt)));
end


toc;
end

%%
load Data_dnaA_titration.mat

rel_dnaA = dnaA_list/dnaA_list(16);
rel_IM = smooth(mean_IM/mean_IM(16),5)';

rel_dnaA_exp = (exp_dnaA/(exp_dnaA(7)))';
rel_IM_exp = (exp_IM/exp_IM(7))';

h=figure;
plot(rel_dnaA,rel_IM,'-','linewidth',2); hold all
plot(exp_dnaA/(exp_dnaA(7)),exp_IM/exp_IM(7),'o','linewidth',2);
xlabel('rel. dnaA expression level');
ylabel('rel. IM');
set(gca,'fontsize',14);
set(gca,'xscale','log');
legend('sim titration-switch','exp data');
xlim([0.1 10]);
