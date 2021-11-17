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
param_stim.F_DBS=[1,2,3,4,5,10,20, 30, 40 ,50,100];
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

%% STEP 3: simulate responses of fit GLM

plotFlag = 1; % plot simulated data
saveFlag = 1; % save simulated data
runs = 1;    % number of trials to simulate

[y, stimcurr, hcurr, r] = simulate_glm(I,dt,k,h,dc,runs,softRect,plotFlag);

figure; plot(y); hold on; plot(stimcurr); plot(hcurr); plot(stimcurr+hcurr)
figure; plot([1:length(I)]*dt,KernelPSTH(y,100,.1,1)); hold on; plot([1:length(I)]*dt,KernelPSTH(spikes,100,.1,1)); plot([1:length(I)]*dt,1e-4*I-.02,'k'); legend('Izhikevic output', 'GLM output','Stimulation'); xlabel('Time (ms)'); title('Firing Rate Validation')
figure; subplot(1,2,1); plot([1:length(k)]*dt,k); title('stimfilter'); subplot(1,2,2); plot([1:length(h)]*dt,h); title('History')   
figure; plot(stimcurr+.2)
hold on; plot(get_stimcurr(1,A,get_stimulations(param_stim),syn_kernel,k, OFFSET))

%% STEP 4: Inferring Steady-state values from high-frequency input data
EOT=5000/dt;
F_DBS=[1,2,3,4,5,10:20:200];
TW=20;
for idx=1:length(F_DBS)
    F_dbs=F_DBS(idx)
    param_stim.F_DBS=[F_dbs];

    I_stim=get_stimulations(param_stim);
    T=length(I_stim)*dt_s;

    syn_kernel=exp(-[0:dt_s:.1]/t_syn).*(1-exp(-[0:dt_s:.1]/2e-3));
    syn_kernel=syn_kernel/max(syn_kernel);
    Temp=conv(I_stim,syn_kernel,'full');
    I=Temp(1:length(I_stim));

    I=A*I'+OFFSET;
    [v, u, spikes, cid] = simulate_izhikevich(cellType,I,dt,jitter,0,0,fid);
    y=spikes(EOT:end);
    x=get_stimcurr(1,A,I_stim(EOT:end),syn_kernel,k,OFFSET);
    x=x(1:length(y));       
    
%     OffsetDc=OFFSET*conv(ones(size(x)),k,'same')+dc;
%     stimcurr=conv(x,k,'full');    
%     stimcurr=stimcurr(1:length(x));
    hcurr=conv(y,(h),'full');
    hcurr=hcurr(1:length(y));
    g=@(a) get_stimcurr(a,A,I_stim(EOT:end)',syn_kernel',k,OFFSET)+hcurr+dc;
    loss=@(lna)  -y'*(g(exp(lna))+log(dt_s))+ dt_s*sum(exp(g(exp(lna)))); %sum((g(a)-g(a_true(idx))).^2);%

    a_estim(idx)=exp(fminsearch(loss,rand));
    a_true(idx)=max(I_stim(EOT:end));
    
    if sum(F_dbs==[50,100,200])
        [GLMout, stimcurr_, hcurr_, r_] = simulate_glm(I,dt,k,h,dc,runs,softRect,0);
        figure 
        plot(dt*[1:length(GLMout)],KernelPSTH(GLMout,TW,dt,1));
        hold on
        plot(dt*[1:length(spikes)],KernelPSTH(spikes,TW,dt,1));
        plot(dt*[1:length(g(a_true(idx)))],exp(-g(a_true(idx))))
        legend('GLM','Izhekevic','\lambda')
        title(['F_{DBS}=',num2str(F_dbs),'Hz'])
    end
    
end
figure; plot(F_DBS,a_true); hold on; plot(F_DBS,a_estim)


function stimcurr=get_stimcurr(a,A,I,syn_kernel,k, OFFSET)
    I=a*A*I/max(I);
    x=conv(I,syn_kernel,'full');
    x=x(1:length(I))+OFFSET;
    stimcurr=conv(x,flip(k),'full');
    stimcurr=stimcurr(1:length(x));
end



