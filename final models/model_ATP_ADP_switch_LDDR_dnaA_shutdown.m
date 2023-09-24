% this script describe the hns dnaA initial model no division

% close all ;
clear all ;

tic;

%% global parameters
DT = 0.5;     %(h)
TT = DT*30 ; % total calculation time (h)
dt = 0.001 ;    % calculation time step (h)
T_shutdown = DT*20 ;

dTdatA = 0 ; % time delay between OriC and datA (h)
dTdars1 = 0.1 ; % time delay between OriC and datA (h)
dTdars2 = 0.2 ; % time delay between OriC and datA (h)
dTdnaA = 0 ; % time delay between OriC and dnaA (min)
CAcc = 0.75 ; % critical cytoplasmic dnaA concentration to initiate DNA replication 

Cperiod = 0.7; 
% repR = 1/Cperiod ; % DNA replication rate = 1 / C period
DNAlength = round(Cperiod/dt) ; % DNA length is set to fix V fork = 1 ;
lambda0 = log(2)/DT ; % cell mass growth rate (h^-1)

D_total = 400;
KDP = 300/D_total;
beta_datA = 600/D_total;
beta_rida = 500/D_total;
alphal = 750/D_total;
alphad1 = 100/D_total;
alphad2 = 643/D_total;

lambda = lambda0;
tseq = 15/60; % hr

%% initiation 
Tlist = (1 : TT/dt)*dt ; % unit [min]
Vlist = nan(1,TT/dt) ;
fATPlist = nan(1,TT/dt) ;
DNAcopy = ones(TT/dt,DNAlength); % starting DNAcopy = 1 ;
DNAcopy(1,1) = 2 ;
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
fATP = 0.8;

% rpt = NaN(1,ceil(TT/DT)+10) ; % DNA replication time
% tert = NaN(1,ceil(TT/DT)+10) ; % DNA termination time
Prep = NaN(1,ceil(TT/DT)+10) ; % position of DNA replication  
Vi =  NaN(1,ceil(TT/DT)+10) ; % initiation mass
Vt =  NaN(1,ceil(TT/DT)+10) ; % termination mass
S = 1;
for t = 1 : Tlist(end)/dt
    
    V = V + lambda*V*dt ;
    DNAlist(t) = sum(DNAcopy(t,:))/DNAlength ;
    No = DNAcopy(t,dTdnaA/dt+1); % copy number of dnaA
    Ndat = DNAcopy(t,dTdatA/dt+1); % copy number of datA
    Nd1 = DNAcopy(t,dTdars1/dt+1); % copy number of dars1
    Nd2 = DNAcopy(t,dTdars1/dt+1); % copy number of dars2

    if t*dt<T_shutdown
        fATP = fATP + dt * ((alphal + alphad1*Nd1/V + alphad2*Nd2/V) * (1-fATP)/(KDP+1-fATP) - (beta_datA + beta_rida) * No/V * fATP/(KDP+fATP) + lambda*(1-fATP));
    else
        fATP = fATP + dt * ((alphal + alphad1*Nd1/V + alphad2*Nd2/V) * (1-fATP)/(KDP+1-fATP) - (beta_datA + beta_rida) * No/V* fATP/(KDP+fATP));
    end
    
    DNAcopy(t+1,:) = DNAcopy(t,:) ;
    % check new replication initiation condition
    if fATP > CAcc && RpO == 1
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

    fATPlist(t) = fATP;
    Vlist(t) = V ;
    OriClist(t) = DNAcopy(t,1) ;
%     rptlist(t,:) = rpt ; 
%     tertlist(t,:) = tert ; 
%     Preplist(t,:) = Prep ; 
%     NAlist(t) = NA; 
%     Ndlist(t) = Ndat; 
    nreplist(t) = nrep;
    RpOlist(t) = RpO; 

end

%% plot
figure('position',[100 100 300 500]);

subplot(2,1,1)
semilogy(Tlist-T_shutdown, Vlist/Vlist(T_shutdown/dt),'linewidth',2); hold all
semilogy(Tlist-T_shutdown, OriClist/OriClist(T_shutdown/dt),'linewidth',2);
plot([0 0], [10^-2 10^2],'k--','linewidth',2);
legend('V','OriC','location','best');
ylabel('Rel. Value');
xlim([-3 3]);
set(gca,'fontsize',14);

subplot(2,1,2)
plot(Tlist-T_shutdown, fATPlist,'linewidth',2); hold all
xlabel('time (hr)');
ylabel('f');
plot([0 0], [0 1],'k--','linewidth',2);
xlim([-3 3]);
ylim([0 1]);
set(gca,'fontsize',14);

