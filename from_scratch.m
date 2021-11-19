    %% Example GLM fitting to Izhikevich neuron
%
% This code simulates data from an Izhikevich neuron, fits a GLM to it,
% simulates responses from that GLM, and then plots a comparison of the 
% simulated GLM responses to the original data.
%
% This code requires the following functions:
%   fminunc (in MATLAB's Optimization Toolbox)
%   generate_izhikevich_stim
%   simulate_izhikevich
%   fit_glm
%   simulate_glm
%   makeBasis_StimKernel
%   makeBasis_PostSpike
%   negloglike_glm_basis
%   negloglike_glm_basis_softrect
%   compare_glm_to_iz
%   normalizecols
%   sameconv
%   logexp1
% All except the first (fminunc) are provided in the package that contains
% this script.


%% STEP 0: Generate DBS inputs
dt_s=1e-4;%sec
% dt=.1; %msec
f=.3;%.1;
U=.1;%0;
F=.03;%.1;
D=.04;%.1;
t_syn=5e-3;
%%%%%
paramTM.f=.3;
paramTM.U=.1;
paramTM.F=.03;
paramTM.D=.04;
paramTM.t_syn=t_syn;
paramTM0=paramTM;
paramTM0.t_syn=1e-6;

T=10;
param_stim.T=T;
param_stim.dt=dt_s;
param_stim.mode='DBS';%% DBS or Poiss or Rand
param_stim.F_DBS=[1,2,3,4,5,10];
param_stim.F_poiss=20;
param_stim.TM=paramTM0;

I_stim=get_stimulations(param_stim);
T=length(I_stim)*dt_s;

syn_kernel=exp(-[0:dt_s:.1]/t_syn).*(1-exp(-[0:dt_s:.1]/2e-3));
syn_kernel=syn_kernel/max(syn_kernel);
Temp=conv(I_stim,syn_kernel,'full');
I_stim=Temp(1:length(I_stim));


%% STEP 1: simulate data from Izhikevich neuron

cellType = 1; % type of Izhikevich neuron (numbered as in Izhikevich 2004)
%   Choose from:
%       1. tonic spiking
%       2. phasic spiking
%       3. tonic bursting
%       4. phasic bursting
%       5. mixed mode
%       6. spike frequency adaptation
%       7. Class 1
%       8. Class 2
%       9. spike latency
%       10. subthreshold oscillations -- not available
%       11. resonator
%       12. integrator
%       13. rebound spike
%       14. rebound burst
%       15. threshold variability
%       16. bistability
%       17. depolarizing after-potential -- not available
%       18. accomodation
%       19. inhibition-induced spiking
%       20. inhibition-induced bursting
%       21. bistability 2 (Not in original Izhikevich paper)

plotFlag = 1; % plot simulated response
saveFlag = 0; % save data to fid, in new folder
fid = pwd;    % root directory for project
T = 10000;    % max time (in ms)
if cellType == 7 || cellType == 8
    T = 20000;  % these behaviors use multiple step heights, so generate more data
end
jitter = 0;   % amount of jitter to add to spike times,
              %   uniformly distributed over [-jitter,jitter], measured in ms



[I, dt] = generate_izhikevich_stim(cellType,T);
if dt~=.1
    error("check the dt")
end
A=200;
OFFSET=5;
I=A*I_stim'+OFFSET;
[v, u, spikes, cid] = simulate_izhikevich(cellType,I,dt,jitter,plotFlag,saveFlag,fid);
subplot(2,1,1); ylim([0 50])
%%            
rate_=KernelPSTH(spikes,100, .1, 1);
figure; plot(rate_)

%% STEP 2: fit GLM

% first choose parameters for basis vectors that characterize the
% stimulus and post-spike filters

% Note that the parameters for basis vectors here were used to fit the
% regular spiking behavior in Weber & Pillow 2017 (Neural Computation).
% For other cell types, different sets of basis vectors were used, so
% simply changing the cell type above will not directly reproduce results
% from the paper.


%%% basis functions for stimulus filter
nkt = 100; % number of ms in stim filter
kbasprs.neye = 0; % number of "identity" basis vectors near time of spike;
kbasprs.ncos = 7; % number of raised-cosine vectors to use
kbasprs.kpeaks = [.1 round(nkt/1.2)];  % position of first and last bump (relative to identity bumps)
kbasprs.b = 10; % how nonlinear to make spacings (larger -> more linear)
%%% basis functions for post-spike kernel
ihbasprs.ncols = 7;  % number of basis vectors for post-spike kernel
ihbasprs.hpeaks = [.1 50];  % peak location for first and last vectors, in ms
ihbasprs.b = 10;  % how nonlinear to make spacings (larger -> more linear)
ihbasprs.absref = 1; % absolute refractory period, in ms

softRect = 0;    % use exponential nonlinearity
plotFlag = 1;    % plot fit
saveFlag = 1;    % save fit to fid, in new folder
maxIter = 300;  % max number of iterations for fitting, also used for maximum number of function evaluations(MaxFunEvals)
tolFun = 1e-12;  % function tolerance for fitting
L2pen = 0;       % penalty on L2-norm of parameter coefficients

[k, h, dc, prs, kbasis, hbasis] = fit_glm(I,spikes,dt,nkt,kbasprs,ihbasprs,[],softRect,plotFlag,maxIter,tolFun,L2pen);

%%
EOT=3000/dt;
F_DBS=100;%[10,20,50:50:200];
for idx=1:length(F_DBS)
F_dbs=F_DBS(idx);
param_stim.F_DBS=[F_dbs];
I_stim=get_stimulations(param_stim);
T=length(I_stim)*dt_s;

syn_kernel=exp(-[0:dt_s:.1]/t_syn).*(1-exp(-[0:dt_s:.1]/2e-3));
syn_kernel=syn_kernel/max(syn_kernel);
Temp=conv(I_stim,syn_kernel,'full');
I_stim=Temp(1:length(I_stim));

TW=500;
I=A*I_stim'+OFFSET;
[v, u, spikes_full, cid] = simulate_izhikevich(cellType,I,dt,jitter,0,0,fid);
figure; plot(KernelPSTH(spikes_full',TW,dt,1)); 

[y, stimcurr, hcurr, r] = simulate_glm(I,dt,k,h,dc,1,softRect,0);
[lambda,r_, g]=GLMval(1*I',y,k,h,dc,dt);
E=@(a) rate_error(a,I,spikes_full,k,h,dc,dt);
rate=KernelPSTH(y,TW,dt,1);
rate_=KernelPSTH(r_,TW,dt,1);
I_=A*I_stim(EOT:end)'+OFFSET;
[v, u, spikes_ss, cid] = simulate_izhikevich(cellType,I_,dt,jitter,0,0,fid);
hold on; plot([EOT:length(spikes_full)],KernelPSTH(spikes_ss',TW,dt,1))
hold on; plot(real(rate_),'K','LineWidth',2); plot(rate); legend('Full Izhekevic','steady-state Izhekevic','full GLMval','full GLMsim')
title(num2str(F_dbs))
% xlim([20000,100000])
% a_estim(idx)=fminsearch(E,rand);
TW = 25;
rate_GLM = KernelPSTH(y,TW,dt,1);
rate_GLM = rate_GLM/mean(rate_GLM) * ave_firing_GLM ;
rate_Izh = KernelPSTH(spikes_full,TW,dt,1);
rate_Izh = rate_Izh/mean(rate_Izh) * ave_firing_Izh ;
figure; hold on,
plot(rate_Izh,'k')
plot(rate_GLM,'r')

end
figure; plot(F_DBS,a_estim)
% figure; plot(rate)
% hold on
% plot(rate_)
% legend('simulated','calculated')

%% helper functions
function [stimcurr,x]=get_stimcurr(a,A,I,EOT,syn_kernel,k, OFFSET)
    I=I(EOT-1e4:end);
    I=a*A*I/max(I);
    x=conv(I,syn_kernel,'full');
    x=[x(1e4+2:length(I));x(end)]+OFFSET;
    stimcurr=conv(x,flip(k),'full');
    stimcurr=stimcurr(1:length(x));
end

function [stimcurr,x]=get_stimcurr2(a,A,I,EOT,syn_kernel,k, OFFSET)
    I=I(EOT-1e4:end);
    I=a*A*I/max(I);
    x=conv(I,syn_kernel,'full');
    x=[x(1e4+2:length(I));x(end)]+OFFSET;
    stimcurr=conv(x,flip(k),'full');
    stimcurr=stimcurr(1:length(x));
end
function [x]=synaptic_transmission(A,I,EOT,syn_kernel)
    I=A*I;
end
function [x]=get_stim_raw(a,A,I,EOT,syn_kernel, OFFSET)
%     I=I(EOT-1e4:end);
    I=a*A*I;
    x=conv(I,syn_kernel,'full');
    x=[x(EOT:length(I));x(end)]+OFFSET;
end
function e=rate_error(a,I,y,k,h,dc,dt)

    [lambda, r, g]=GLMval(a*I,y,k,h,dc,dt);
    TW=10;
    rate=KernelPSTH(y,TW,dt,1);
    rate_=KernelPSTH(r,TW,dt,1);
    e=mean((rate-rate_).^2);
    
end
function [lambda, rate, g]=GLMval(x,y,k,h,dc,dt)
    refreshRate = 1000/dt; 
   stimcurr=sameconv(x,k);
   hcurr   =sameconv(y,flip(h));
   g=stimcurr+hcurr+dc;
   lambda=exp(g);%1-exp(-g/refreshRate);
   rate=1-exp(-lambda*(dt*1e-3));
end
function y=CUT(x,END)
    y=x(1:END);
end
