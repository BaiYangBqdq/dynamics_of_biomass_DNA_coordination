% this script describe the hns dnaA initial model + adder division model
% the core of this model is OriC/determine adder at each division

% close all ; 
clear all ;


tic;
dnaA_list = 2.^[0:0.2:6];
mean_G_conc = nan(size(dnaA_list));
for id = 1:length(dnaA_list)
%% global parameters
TT = 1000 ; % total calculation time (min)
dt = 0.1 ; % calculation time step (min)
dtR = 1 ; % time step for recording (min)
RTlist = dtR : dtR : TT ; % time list (min)
% Ctime = 250; % changing time of hns expression (min)


CAcc = 10 ; % critical cytoplasmic dnaA concentration to initiate DNA replication 
alpha0 = dnaA_list(id);
datAn = 100 ; % dnaA binding site on datA
Chns = 300 ; % number of hns binding to each chromosome
alphaH0 = 0;
Tseq = 10;
dTdatA = 0 ; % time delay between OriC and datA (min)

DT = 27; 
Cperiod = 42;
Tcycle = Cperiod*3/2;
DNAlength = ceil(Cperiod/dt) ; % DNA length is set to fix V fork = 1 ;
lambda0 = log(2)/DT ; % cell mass growth rate (min^-1)


CNi = 1000; % initial cell number
CNmax = 10000; % maxi cell number
CNdilut = 0.8*CNmax ; % dilute the cells at this cell number
dilutRatio = 10 ; 
Tdilut = nan(1,20) ; % diluting ratio % dilution time

Tlist = dt:dt:TT;
alpha0_list = alpha0*ones(size(Tlist));


%% single cell phisiological properties structure

A = nan ; % cellular dnaA protein ; 
Af = nan ; % cytoplasmic dnaA protein
H = nan ; % cytoplasmic hns protein
Abox = nan;
V = nan ; % cell volume ;

DNAcopy = ones(1,DNAlength); % DNA copy;
OriC = DNAcopy(1); % cellular OriC
Nd = DNAcopy(dTdatA/dt+1); % copy number of datA
DNA = sum(DNAcopy)/DNAlength ; % cellular DNA amount

nrep = nan ; % number of DNA replication. starting from 1, not 0
ndiv = nan;
lambda = nan; 
alpha = nan; 

rpt = NaN(1,ceil(TT/DT)+10) ;
divt =  NaN(1,ceil(TT/DT)+10)  ; % division time
Prep = NaN(1,ceil(TT/DT)+10)  ; % position of DNA replication  
Tdiv = [inf];

% Cell = repmat(C,CNmax, 1) ;
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
% divt =  NaN(1,ceil(TT/DT)+10)  ; % division time
% Prep = NaN(1,ceil(TT/DT)+10)  ; % position of DNA replication  
% Tdiv = [inf];


%% main loop
CN = 1 ;
ndilut = 0 ; % number of dilution times

for t = 1 : length(Tlist)
    
    alpha0 = alpha0_list(t);
    
    OriC = DNAcopy(1);
    DNA = sum(DNAcopy)/DNAlength ;
    Nd = DNAcopy(dTdatA/dt+1); % copy number of datA
    A = A + alpha*V*dt;
    V = V + lambda*V*dt ;
    H = alphaH0*V;
    Abox = max((Nd*datAn + DNA*Chns - H),0);
    Af = max(A-Abox,0);
    

    
    % check new replication initiation condition
%     if Af/V > CAcc && RpO == 1
    if Af/V > CAcc
        nrep = nrep +1; % replication number +1
        Prep(nrep) = 1 ; % generate new replication position
        rpt(nrep) = t*dt; % replication time +dt
        RpO = 0; % stop new replicaiton initiation
    end
    
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

rpt(isnan(rpt))=[];
G_conc = Rec_DNA./Rec_V;
IM = Rec_V./Rec_OriC;
mean_G_conc(id) = mean(G_conc(rpt(end-1)/dtR:rpt(end)/dtR));
mean_IM(id) = mean(IM(rpt(end-1)/dtR:rpt(end)/dtR));
end

toc;
%% plot time course data
% 
% ht=figure(); 
% semilogy(RTlist(1:600)/60,Rec_V(1:600),'linewidth',2); hold all
% semilogy(RTlist(1:600)/60,Rec_OriC(1:600),'linewidth',2); hold all
% % semilogy(RTlist/60,Rec_DNA,'linewidth',2); hold all
% % semilogy(RTlist/60,Rec_A,'linewidth',2); hold all
% xlabel('Time (hr)');
% legend('mass','oriC','DNA','dnaA','location','NorthWest'); hold off
% set(gca,'fontsize',14);
% drawnow
% 
% h=figure(); 
% subplot(2,1,1)
% plot(Tlist(3000:6000)/60,list_OriC(3000:6000)./list_V(3000:6000),'linewidth',2); hold all
% ylabel('oriC / mass');set(gca,'fontsize',14);
% subplot(2,1,2)
% plot(Tlist(3000:6000)/60,list_Af(3000:6000)./list_V(3000:6000),'linewidth',2); hold all
% xlabel('Time (hr)');
% ylabel('free dnaA / mass');set(gca,'fontsize',14);
% drawnow

% figure;
% plot(dnaA_list,mean_G_conc,'o-','linewidth',2)
% xlabel('dnaA expression level');
% ylabel('Chromosome concentration');
% set(gca,'fontsize',14);


h=figure;
plot(dnaA_list/dnaA_list(16),smooth(mean_IM/mean_IM(16),5),'-','linewidth',2); hold all
load Data_dnaA_titration.mat;
plot(exp_dnaA/(exp_dnaA(7)),exp_IM/exp_IM(7),'o','linewidth',2);
xlabel('rel. dnaA expression level');
ylabel('rel. IM');
set(gca,'fontsize',14);
set(gca,'xscale','log');
legend('sim with hns','exp data');

% saveas(h,'IM_hns_fit_hns2000','fig');

%% save data
rel_dnaA = dnaA_list/dnaA_list(16);
rel_IM = smooth(mean_IM/mean_IM(16),5)';

rel_dnaA_exp = (exp_dnaA/(exp_dnaA(7)))';
rel_IM_exp = (exp_IM/exp_IM(7))';

% save('IM_hns_fit_hns0.mat','rel_dnaA','rel_IM','rel_dnaA_exp','rel_IM_exp');
