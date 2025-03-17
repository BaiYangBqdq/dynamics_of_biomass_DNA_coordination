% this script describe the hns dnaA initial model + adder division model
% the core of this model is OriC/determine adder at each division

close all ; clear all ;

hns_list = 50*2.^[0:0.2:6];
mean_IM = nan(size(hns_list));
mean_DnaA_activity = nan(size(hns_list));

tic;
for id = 1:length(hns_list)
%% global parameters
TT = 10000 ; % total calculation time (min)
dt = 0.1 ; % calculation time step (min)
dtR = 1 ; % time step for recording (min)
RTlist = dtR : dtR : TT ; % time list (min)
% Ctime = 500; % changing time of hns expression (min)


CAcc = 10 ; % critical cytoplasmic dnaA concentration to initiate DNA replication 
alpha0 = 5;
datAn = 100 ; % dnaA binding site on datA
Chns = 300 ; % number of hns binding to each chromosome
alphaH0 = hns_list(id);
Tseq = 10;
dTdatA = 10 ; % time delay between OriC and datA (min)

DT = 27; 
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
Vi = NaN(1,ceil(TT/DT)+10) ;
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
        Vi(nrep) = V/DNAcopy(1);
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




Vi(isnan(rpt))=[];
rpt(isnan(rpt))=[];
mean_IM(id) = mean(Vi(10:end));
mean_DnaA_activity(id) = mean(list_Af(round(rpt(100:end)/dt))./list_V(round(rpt(100:end)/dt)));


%% plot time course data
% 
% ht=figure(); 
% semilogy(RTlist/60,Rec_V/(Rec_V(Ctime/dtR)),'linewidth',2); hold all
% semilogy(RTlist/60,Rec_OriC/(Rec_OriC(Ctime/dtR)),'linewidth',2); hold all
% semilogy([Ctime Ctime]/60,[10^-10 10^5],'k--','linewidth',2); hold all
% xlim([0 15]); ylim([10^-10 10^5])
% xlabel('Time (hr)'); ylabel('Rel. value');
% legend('mass','oriC','location','NorthWest'); hold off
% set(gca,'fontsize',14);
% drawnow
% 
% 
%%
h=figure('position',[100 100 500 600]); 
% subplot(3,1,1)
% plot((Tlist(3000:7000)-Ctime)/60,alpha0_list(3000:7000),'linewidth',2); hold all
% ylabel('\alpha_A'); xlim([-3 3]);
% set(gca,'fontsize',14);
subplot(2,1,1)
plot(Tlist(round(rpt(100)/dt):round(rpt(end)/dt))/60,list_Af(round(rpt(100)/dt):round(rpt(end)/dt))./list_V(round(rpt(100)/dt):round(rpt(end)/dt)),'linewidth',2); hold all
plot(Tlist(round(rpt(100)/dt):round(rpt(end)/dt))/60, 10*ones(size((Tlist(round(rpt(100)/dt):round(rpt(end)/dt))/60))),'k--','linewidth',2);
ylabel('free dnaA / mass');
set(gca,'fontsize',14);
subplot(2,1,2)
plot(Tlist(round(rpt(100)/dt):round(rpt(end)/dt))/60,list_OriC(round(rpt(100)/dt):round(rpt(end)/dt))./list_V(round(rpt(100)/dt):round(rpt(end)/dt)),'linewidth',2); hold all
ylabel('oriC / mass');
% xlim([0 6]);

xlabel('Time (hr)');
set(gca,'fontsize',14);
drawnow


end
%%
rel_hns = hns_list/hns_list(16);
rel_IM = mean_IM/mean_IM(16);
rel_DnaA = mean_DnaA_activity/mean_DnaA_activity(16)';

h=figure;
plot(rel_hns,rel_IM,'-','linewidth',2); hold all
plot(rel_hns,rel_DnaA,'-','linewidth',2); hold all
xlabel('rel. hns expression level');
ylabel('rel. value');
set(gca,'fontsize',14);
set(gca,'xscale','log');
legend('IM','DnaA activity');
xlim([0.1 10]);