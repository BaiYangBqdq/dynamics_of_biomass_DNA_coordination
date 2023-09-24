% this script describe the hns dnaA initial model no division

close all ;
clear all ;

%% global parameters
TT = 70000 ; % total calculation time (min)
dt = 0.2 ; % calculation time step (min)
% Tlist = dt:dt:TT;
% Aflist = nan(size(Tlist));

CAcc = 10 ; % critical cytoplasmic dnaA concentration to initiate DNA replication 
alpha0 = 1;
datAn = 100 ; % dnaA binding site on datA
Chns = 300 ; % number of hns binding to each chromosome
alphaH0 = 0;
Tseq = 0;
dTdatA = 0 ; % time delay between OriC and datA (min)

DT = 60; 
Cperiod = 42;
DNAlength = ceil(Cperiod/dt) ; % DNA length is set to fix V fork = 1 ;
lambda0 = log(2)/DT ; % cell mass growth rate (min^-1)

epsilon_A = 0.1;
epsilon_V = 0.1;
Prep = 1;

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
Oi = NaN(1,ceil(TT/DT)+10) ;
%% main loop
for t = 1 : TT/dt

    
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
        Vi(nrep) = V; % replication time +dt
        Oi(nrep) = OriC;
%         RpO = 0; % stop new replicaiton initiation
        alpha = alpha0*(epsilon_A * randn(1)+1);
        lambda = lambda0*(epsilon_V * randn(1)+1);
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
    
    
%     Aflist(t) = Af/V;
end



%% adder plot
dbin = 0.005;
Vt = rpt;
Vip = Vi(10:end-1)./Oi(10:end-1);
dVip = Vi(11:end)./Oi(11:end)*2-Vi(10:end-1)./Oi(10:end-1);
dVit = rpt(11:end)-rpt(10:end-1);
Vip(isnan(dVip))=[];
dVit(isnan(dVip))=[];
dVip(isnan(dVip))=[];


h5 = figure('position',[100 400 900 300]); 
subplot(1,3,2)
bin = min(Vip):dbin:max(Vip);
disc = discretize(Vip,bin);
for i = 1:length(bin)
    meanVip(i) = nanmean(dVip(disc==i));
    stddVip(i) = nanstd(dVip(disc==i));
end
scatter(Vip,dVip,10,'filled'); hold all
% errorbar(bin,meanVip,stddVip,'ro-','linewidth',2);
% xlim([min(Vip)*0.99 max(Vip)*1.01]);
% ylim([min(Vip)*0.99 max(Vip)*1.01]);
pa=polyfit(Vip,dVip,1);
sla = pa(1);
xlabel('initiation mass');
ylabel('added mass between initiations');
title(sprintf('slope = %.2f',sla));
set(gca,'fontsize',14);

subplot(1,3,1)
Vip1 = Vip(1:end-1);
Vip2 = Vip(2:end);
bin = min(Vip1):dbin:max(Vip1);
disc = discretize(Vip1,bin);
for i = 1:length(bin)
    meanVip2(i) = nanmean(Vip2(disc==i));
    stdVip2(i) = nanstd(Vip2(disc==i));
end
scatter(Vip1,Vip2,10,'filled'); hold all
% errorbar(bin,meanVip2,stdVip2,'ro-','linewidth',2);
% xlim([min(Vip)*0.99 max(Vip)*1.01]);
% ylim([min(Vip)*0.99 max(Vip)*1.01]);
pb=polyfit(Vip(1:end-1),Vip(2:end),1);
slb = pb(1);
xlabel('initiation mass');
ylabel('next initiation mass');
title(sprintf('slope = %.2f',slb));
set(gca,'fontsize',14);

subplot(1,3,3)
bin = min(Vip):dbin:max(Vip);
disc = discretize(Vip,bin);
for i = 1:length(bin)
    meandVit(i) = nanmean(dVit(disc==i));
    stddVit(i) = nanstd(dVit(disc==i));
end
scatter(Vip,dVit,10,'filled'); hold all
% errorbar(bin,meandVit,stddVit,'ro-','linewidth',2);
% xlim([min(Vip)*0.99 max(Vip)*1.01]);
% ylim([min(dVit)*0.99 max(dVit)*1.01]);
p=polyfit(Vip(1:end-1),dVit(2:end),1);
slc = p(1);
xlabel('initiation mass');
ylabel('time between initiation ');
title(sprintf('slope = %.2f',slc));
set(gca,'fontsize',14);


%% adder plot normalized

meanVip=nanmean(Vip);       scaleVip = Vip/meanVip;
meandVit=nanmean(dVit);     scaledVit = dVit/meandVit;
meandVip=nanmean(dVip);     scaledVip = dVip/meandVip;
bin = 0.5:0.05:1.5;

h5 = figure('position',[100 400 900 300]); 
subplot(1,3,2)
disc = discretize(scaleVip,bin);
for i = 1:length(bin)
    meanscaleVip(i) = nanmean(scaledVip(disc==i));
    stddscaleVip(i) = nanstd(scaledVip(disc==i));
end
scatter(scaleVip,scaledVip,10,'filled'); hold all
plot(bin,meanscaleVip,'ro-','linewidth',2);
xlim([0.5 1.5]);    ylim([0.5 1.5]);
pa=polyfit(scaleVip,scaledVip,1);
sla = pa(1);
xlabel('scaled init. mass');
ylabel('scaled added init. mass');
title(sprintf('slope = %.2f',sla));
set(gca,'fontsize',14);

subplot(1,3,1)
scaleVip1 = scaleVip(1:end-1);
scaleVip2 = scaleVip(2:end);
disc = discretize(scaleVip1,bin);
for i = 1:length(bin)
    meanscaleVip2(i) = nanmean(scaleVip2(disc==i));
    stdscaleVip2(i) = nanstd(scaleVip2(disc==i));
end
scatter(scaleVip1,scaleVip2,10,'filled'); hold all
plot(bin,meanscaleVip2,'ro-','linewidth',2);
xlim([0.5 1.5]);
ylim([0.5 1.5]);
pb=polyfit(scaleVip1,scaleVip2,1);
slb = pb(1);
xlabel('scaled init. mass');
ylabel('scaled next init. mass');
title(sprintf('slope = %.2f',slb));
set(gca,'fontsize',14);

subplot(1,3,3)
disc = discretize(scaleVip,bin);
for i = 1:length(bin)
    meanscaledVit(i) = nanmean(scaledVit(disc==i));
    stdscaledVit(i) = nanstd(scaledVit(disc==i));
end
scatter(scaleVip,scaledVit,10,'filled'); hold all
plot(bin,meanscaledVit,'ro-','linewidth',2);
xlim([0.5 1.5]);
ylim([0.5 1.5]);
p=polyfit(scaleVip,scaledVit,1);
slc = p(1);
xlabel('scaled init. mass');
ylabel('scaled inter init. time');
title(sprintf('slope = %.2f',slc));
set(gca,'fontsize',14);
