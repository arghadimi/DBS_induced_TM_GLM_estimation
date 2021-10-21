%% Simulation parameters
close all
clear all
clc

dt=1e-4;

f=.3;%.1;
U=.1;%0;
F=.03;%.1;
D=.04;%.1;
t_syn=3e-3;
%%%%%
paramTM.f=f;
paramTM.U=U;
paramTM.F=F;
paramTM.D=D;
paramTM.t_syn=t_syn;
paramTM0=paramTM;
paramTM0.t_syn=1e-6;

I_inf=[];
I_inf_comp=[];
for i=10:1:200
    I_inf(end+1)=SteadyStateCurrent(paramTM,i)./DTM_DBS(paramTM,i,1);
    temp=DTM_DBS(paramTM,i,100)./DTM_DBS(paramTM,i,1);
    I_inf_comp(end+1)=temp(100);
end

figure; plot([DTM_DBS(paramTM,5,60)]); hold on ;plot([DTM_DBS(paramTM,100,60)]);
legend('1Hz', '100 Hz')
title("sample STP response")
%% DESINGING the basis functions

figure; ax=subplot(1,1,1); title('Designed Coupling Filter')
coupling.delay=.1/dt;
coupling.fie=pi*[0:.5:3];
coupling.dt=dt;
coupling.c=1;
coupling.a=200;
coupling.basis=RCosBasis(coupling);
% coupling.bta=[.7,.8,1.3,1.5,.9,.7,.5]/15;
coupling.bta=[.3,.4,.7,.9,.3,.15,.1]/20;% if there is no history its better to use a filter like this
coupling.filter=coupling.bta*coupling.basis;
coupling.nfilt=size(coupling.basis,1);
hold on
plot(1e3*dt*[1:coupling.delay],coupling.bta.*coupling.basis','r');
plot(1e3*dt*[1:coupling.delay],(coupling.filter),'k','LineWidth',2);
set(ax,'XTick',[0, 50 , 100])
xlabel Time(ms)
title 'Designed History Filter'
%


hist.delay=0;%.02/dt;
hist.fie=pi*[1,1.5];
hist.dt=dt;
hist.c=1;
hist.a=2000;
hist.basis=RCosBasis(hist);

if hist.delay>0
    hist.bta=[-5,-1];
    hist.filter=hist.bta*hist.basis;
    hist.nfilt=size(hist.basis,1);
else
    hist.bta=[];
    hist.filter=0;%hist.bta*hist.basis;
    hist.nfilt=size(hist.basis,1);
end

if hist.delay>0
    figure; ax=subplot(1,1,1); title('Designed History Filter')
    hold on
    plot(1e3*dt*[1:hist.delay],hist.bta.*hist.basis','r');
    plot(1e3*dt*[1:hist.delay],(hist.filter),'k','LineWidth',2);
    set(ax,'XTick',[0, 10 , 20])
    xlabel Time(ms)
    title History
end
OFFSET=2;
%% getting the synthetic input patterns
T=10;
param_stim.T=T;
param_stim.dt=dt;
param_stim.mode='DBS';
param_stim.F_DBS=[1,2,3,4,5,10];
param_stim.F_poiss=20;
param_stim.TM=paramTM0;

I_stim=get_stimulations(param_stim);
T=length(I_stim)*dt;

syn_kernel=exp(-[0:dt:.1]/t_syn);
syn_kernel=syn_kernel/max(syn_kernel);
Temp=conv(I_stim,syn_kernel,'full');
I_stim=Temp(1:length(I_stim));

figure
plot(dt*[1:length(I_stim)],I_stim)
title PSC_{stim}
xlabel Time(sec)

k=conv(I_stim,coupling.filter,'full')+OFFSET;
V_bin=zeros(1,hist.delay);


%% Generating a Poisson firing train
for i=hist.delay+1:length(I_stim)
    if hist.delay>0
        h(i)=(V_bin(i-floor(hist.delay):i-1))*flip(hist.filter');
    else
        h(i)=0;
    end
    y_rate(i)=exp(k(i)+h(i));
    V_bin(i)=rand<1-exp(-y_rate(i)*dt) && sum(V_bin(i-5:i-1))==0;
end
% V_bin=spike_poiss2(T,dt,y_rate);
p=[];
figure; p(1)=ScatterPlot(V_bin,dt,-.5,-1,'k'); hold on
p(2)=plot(dt:dt:T,I_stim,'r','LineWidth',2);
xlim([0 10])

xlabel('Time(sec)')
legend(p,{'Output Spikes','Input Stimulations'})


%% Constructing the design Matrix

XY=zeros(length(I_stim),coupling.nfilt);
for i=1:coupling.nfilt
    temp=conv([0,I_stim],coupling.basis(i,:));
    XY(:,i)=temp(1:length(XY));
end
V_bin_pad=[zeros(1,99),V_bin(1:end-99)];
XX=zeros(length(V_bin_pad),hist.nfilt); %% There is the problem that must be double-checke
for i=1:hist.nfilt
    temp=conv(V_bin,hist.basis(i,:),'full');
    XX(:,i)=temp(1:length(V_bin_pad));
end
%% GlmFit 
if hist.delay>0
    Xdsgn =[XY,XX];
else 
   Xdsgn =[XY];
end
y=conv(V_bin,hann(40),'full'); % smoothing
y=y(1:length(V_bin));

AvgFR=sum(V_bin)/T;
wk=3*dt;
rate_ = KernelPSTH(V_bin,wk,dt,1);
rate_=rate_/1330;
% rate_ = rate_*AvgFR/mean((rate_));





bta_estim = glmfit(Xdsgn,[min(rate_,1);ones(size(rate_))]','binomial','link','logit');%, 'constant','off');%,'offset',OFFSET*ones(size(rate_)));
RrME_GLM= glmval(bta_estim,Xdsgn,'logit');%,'constant','off','offset',-.5*ones(size(rate_')));

% bta_estim = glmfit(Xdsgn,[rate_]','poisson','link','log','constant','off','offset',OFFSET*ones(size(rate_)));
% RrME_GLM= glmval(bta_estim,Xdsgn,'log','constant','off','offset',OFFSET*ones(size(rate_')));

coupling.bta_estim=bta_estim(2:coupling.nfilt+1)';
coupling.filter_estim=coupling.bta_estim*coupling.basis;
hist.bta_estim=bta_estim(end-hist.nfilt+1:end)';
hist.filter_estim=hist.bta_estim*hist.basis;

figure
hold on
plot(rate_);
plot(RrME_GLM);
plot(.5*V_bin-1,'k')
hold off
legend('Reference FR','GLM out','Reference binary')

figure
subplot(2,1,1)
plot((coupling.filter_estim))
hold on
plot((coupling.filter))
legend('estimated','original')
subplot(2,1,2)
plot((hist.filter_estim))
hold on
plot((hist.filter))

%% Steady_state estimation
% unibserved frequency
F_DBS=[5,10,20,30,50:50:200];%[5,10,20,30,50,100,200];

clear A_estim XY XX
K=1;
F_idx=0;
V_bin_new=zeros(8,500000);
for F_dbs=F_DBS
% the PSC amplitude that it is to be estimated later
F_idx=F_idx+1;

T_dbs=1/(dt*F_dbs);
t_stim=floor(T_dbs:T_dbs:T/dt);
I_new=zeros(1,round(T/dt));%filter(1,[1 -0.9],I_rand);
I_new(t_stim)=DTM_DBS(paramTM0,F_dbs,length(t_stim))/DTM_DBS(paramTM0,F_dbs,1);
I_SS(F_idx)=I_new(t_stim(end));
Temp=conv(I_new,syn_kernel,'full');
I_new=Temp(1:length(I_new));


k=conv(I_new,coupling.filter,'full')+ OFFSET;
%     V_bin_new(m,:)=zeros(1,hist.delay);

for i=hist.delay+1:length(I_new)
    if hist.delay>0
        h(i)=(V_bin_new(F_idx,i-hist.delay:i-1))*flip(hist.filter');
    else
        h(i)=0;
    end
    y_rate(i)=exp(k(i)+h(i));
    V_bin_new(F_idx,i)=rand<1-exp(-y_rate(i)*dt) && sum(V_bin_new(F_idx,i-5:i-1))==0;
end


% I_new=I_new(20000:end);
% V_bin_new=V_bin_new(20000:end);
EOT=50000;
I_new_cut(F_idx,:)=I_new((EOT):end);
I_new_cut(F_idx,:)=zeros(size(I_new_cut(F_idx,:)));
I_new_cut(F_idx,floor(T_dbs:T_dbs:length(I_new_cut(F_idx,:))))=1;
Temp=conv(I_new_cut(F_idx,:),syn_kernel,'full');
I_new_cut(F_idx,:)=Temp(1:length(I_new_cut(F_idx,:)));

V_bin_new_cut(F_idx,:)=V_bin_new(F_idx,(EOT):end);

XY=zeros(length(I_new_cut(F_idx,:)),coupling.nfilt);
for i=1:coupling.nfilt
    temp=conv([0,I_new_cut(F_idx,:)],coupling.basis(i,:));
    XY(:,i)=temp(1:length(XY));
end

if hist.delay>0
    V_bin_new_pad=[zeros(1,hist.delay-1),V_bin_new_cut(F_idx,1:end-hist.delay+1)];
    XX=zeros(length(V_bin_new_pad),hist.nfilt); %% There is the problem that must be double-checke
    for i=1:hist.nfilt
        temp=conv(V_bin_new_cut(F_idx,:),hist.basis(i,:),'full');
        XX(:,i)=temp(1:length(V_bin_new_pad));
    end

    Xdsgn = [XY,XX];
else 
    Xdsgn=[XY];
end
%
AvgFR=sum(V_bin_new_cut(F_idx,:))/T;
rate_new = KernelPSTH(V_bin_new_cut(F_idx,:),wk,dt,1);
% rate_new = rate_new*AvgFR/mean(rate_new);
rate_new=rate_new/1330;
%bta_estim(1); 
if hist.delay>0
    y_estim=@(A) glmval([bta_estim(1);A*coupling.bta_estim';hist.bta_estim'],Xdsgn,'logit');%,'constant','off','offset',OFFSET);
else
    y_estim=@(A) glmval([bta_estim(1);A*coupling.bta_estim'],Xdsgn,'logit');%,'constant','off','offset',OFFSET);
end
error= @ (lnA) (mean((rate_new'-y_estim(1.*exp(lnA))).^2));
% figure
% plot(rate_new'); hold on; plot(y_estim(I_SS(F_idx))); legend('REF','Estim')
% xlim([0, 1e4])
[LnA_estim,fval(F_idx)]=fminsearch(error,-rand);
%
A_estim(F_idx)=exp(LnA_estim);
end


figure; plot(F_DBS,A_estim);
hold on
plot(F_DBS,I_SS)
clear I_inf
for i=1:length(F_DBS); I_inf(i)=SteadyStateCurrent(paramTM0,F_DBS(i))/DTM_DBS(paramTM,F_DBS(i),1); end
plot(F_DBS,I_inf,'k--')

%% Estimating the TM parameters
ub=[0 0 1 1];
lb=[-5 -5 -5 -5];
x0=rand(1,4).*(ub-lb)+lb;
[fitresult, gof] = Gradient_fit(F_DBS,(DTM_DBS(paramTM,1,1))*A_estim,x0,1e-5 ,ub, lb,500,'AllLog',0)
f_estim=exp(fitresult.f);
U_estim=exp(fitresult.U);
F_estim=exp(fitresult.F);
D_estim=exp(fitresult.D);
%probably need to double check the formula
paramTM_estim.f=f_estim;
paramTM_estim.U=U_estim;
paramTM_estim.F=F_estim;
paramTM_estim.D=D_estim;
paramTM_estim.t_syn=t_syn;
paramTM0_estim=paramTM_estim;
paramTM0_estim.t_syn=1e-6;

for i=1:length(F_DBS); I_inf_estim(i)=SteadyStateCurrent(paramTM0_estim,F_DBS(i))/DTM_DBS(paramTM_estim,F_DBS(i),1); end
plot(F_DBS,I_inf_estim)

legend('GLM infered Steady-state','Correct steady-state','Analythical solution','TM estimated')
%%
yes='y';
no='n';
response=input('Do you want to plot the Mult-frequency time traces?(yes/no)');
if strcmpi(response,yes) || strcmpi(response,'yes') || strcmpi(response,'y')
    figure; 
    for F_idx=1:length(F_DBS)
        F_dbs=F_DBS(F_idx);
        T_dbs=1/(dt*F_dbs);
        t_stim=floor(T_dbs:T_dbs:T/dt);
        I_new=zeros(1,round(T/dt));%filter(1,[1 -0.9],I_rand);
        I_new(t_stim)=DTM_DBS(paramTM0_estim,F_dbs,length(t_stim))/DTM_DBS(paramTM0_estim,F_dbs,1);
        I_SS(F_idx)=I_new(t_stim(end));
        Temp=conv(I_new,syn_kernel,'full');
        I_new=Temp(1:length(I_new));


        k=conv(I_new,coupling.filter,'full')+ OFFSET;
        %     V_bin_new(m,:)=zeros(1,hist.delay);

        for i=hist.delay+1:length(I_new)
            h(i)=(V_bin_new(F_idx,i-hist.delay:i-1))*flip(hist.filter');
            y_rate(i)=exp(k(i)+h(i));
            V_bin_new(F_idx,i)=rand<1-exp(-y_rate(i)*dt) && sum(V_bin_new(F_idx,i-5:i-1))==0;
        end

        rate_new1 = KernelPSTH(V_bin_new_cut(F_idx,:),100*wk,dt,1);
        rate_new2 = KernelPSTH(V_bin_new(F_idx,:),100*wk,dt,1);

        subplot(2,4,F_idx)
        plot(rate_new1); hold on; plot(rate_new2)
        xlim([0, 5e5])
        title([num2str(F_dbs),'Hz'])
    end
end