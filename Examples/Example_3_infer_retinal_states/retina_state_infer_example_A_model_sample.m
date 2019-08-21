%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up workspace
close all; clear all;
addpath('./','../','../../','../../NFCP'); 
NFCP_init

% Set this to the path of NFCP_datasets on your system
% (Contact mrule7404@gmail.com for a copy of NFCP_datasets)
% Here is a dropbox share link (current as of 2019/08/18)
% https://www.dropbox.com/sh/apifkoowouvq2jc/AABLysnAJWs1YY2fmB4HKVRWa?dl=0
datadir  = '/home/mer49/Dropbox (Cambridge University)/NFCP_datasets/';

% Script options
savefigs  = false;
showevery = 10; 
df        = 'P11_binned_100ms_20150219_leptin_retina1_control';
SRES      = 10; % spatial resolution (10, 15 and 20 are available)
Nsample   = 5000;

% Locate data file
subdir   = 'binned_spikes/spatiotemporal_binned/renewal_model/version_4_selected_datasets/';
fn       = sprintf('%d/%s_resolution_%d',SRES,df,SRES);
suffix   = '_preinitialized_using_renewal_model.mat';
filepath = [datadir subdir fn suffix];
loadarchive 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model configuration

% Heuristic model configuration tuned for P11_binned_100ms_20150219_leptin_retina1_control
% One change for sampling: setting a nonzero rate of spontaneous events
model.rAA           = 10;    % Excito-Excitatory interaction strength
model.thr           = 0.10;  % Activation threshold
model.sigma         = 0.1;   % Standard-deviation for excitatory interaction kernel
model.description = [ ...
%    Q  A  R   rate
    -1  1  0   10.0  % spontaneous excitation 
     0 -1  1   1.8   % enter slow refractory loop
     1  0 -1   0.1   % slow refractory recovery
     ];

%{
% Incorrect parameters, what happens to the dynamics? 
model.rAA           = 5;    % Excito-Excitatory interaction strength
model.thr           = 0.1;  % Activation threshold
model.sigma         = 0.1;   % Standard-deviation for excitatory interaction kernel
model.description = [ ...
%    Q  A  R   rate
    -1  1  0   10.0  % spontaneous excitation 
     0 -1  1   1.8   % enter slow refractory loop
     1  0 -1   0.025   % slow refractory recovery
     ];
%}

% Spiking observations are turned off for free-run sampling example
model.gain          = 0;     % Gain for linear Cox-process observation model for spikes
model.bias          = 0;     % Bias for linear Cox-process observation model for spikes
model.alpha         = 1;     % Dispersion paramter, 1=Poisson
model.inhomogeneous_gain = 1;

% Algorithim configuration parameters
model.dt            = 0.1;  % time interval between observations
model.n             = SRES; % Simulation grid size
model.update        = 'Laplace-subspace'; % Measurement update method
model.minrate       = 1e-6;  % Regularization parameter to avoid zeros
model.cutoff        = false; % Remove scales finer than interaction radius?
model.reg_state_var = 1e-4;  % Added uncertainty in state 
model.reg_count_var = 0;     % Added uncertainty in neuron count
model.reg_diag      = 1e-6;  % Added diagonal regularizaer
model.reg_inverse   = 1e-6;  % Diagonal regularizer for matrix inversions
model.ini_state_var = 1e-0;  % Initial variance in state estimate
model.ini_count_var = 0;     % Initial variance for neuron count
model.ini_reg_diag  = 1e-6;  % Initial variance added to diagonal
model.reg_precision = 1e-6;  % Regularization for precision matrix
model.maxiter       = 25;    % Limit Newton-raphson measurement iterations
model.verbosity     = 0;     % How many details to print
model.safety        = 0;     % Whether to do extra numeric stability checks
effectivepopsize    = 100;   % Effective population size for scaling noise

% Colors and rescaling for plotting
model.cscale        = [1 5 1];
model.colors        = [0 0 1; % Q: blue
                       1 0 0; % A: red
                       0 1 0];% R: green
model.names = strsplit('Q A R');
fprintf('model.ss_rescale = %f\n',model.ss_rescale);
ratescale = model.cscale(2);
model.names{2} = sprintf('A×%d',ratescale);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set system size rescaling based on effective population size
model.ss_rescale = 1/effectivepopsize;
model = initializeModel(model);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Soft boundary around active regions
% (messy code gets reflected boundary conditions for a spatial blur)counts = {};
for i=1:size(xydata,2), counts{i} = binCounts(model,xydata{i}); end
counts = cell2mat(counts);
for i=1:N*N, m2(i) = mean(counts(i,1:end)); end
U = m2>prctile(m2,15);
U = reshape(U,N,N);
Upadded = [U fliplr(U)];
Upadded = [Upadded; flipud(Upadded)];
K       = floor(N/2);
Upadded = circshift(Upadded,[K,K]);
Upadded = gaussian2DblurOperator(N*2,0.1*N)*Upadded(:);
Upadded = reshape(Upadded,N*2,N*2);
U       = Upadded(K:K+N-1,K:K+N-1);
U       = U(:)';
ini = [U*0.7 U*0 U*0.3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample from model
fprintf(1,'Sampling from model...\n');
[ini,xydata,rates,simulatedM] = stateSample(model,{...
    'doplot'          ,true 
    'Nsample'         ,Nsample 
    'Nburn'           ,250  
    'upscale'         ,8    
    'skipsim'         ,showevery
    'oversample'      ,10   
    'effectivepopsize',1
    'save_figure'     ,savefigs
    'efraction'       ,0.04
    'ini'             ,ini
    'rescale'         ,true
    });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistical summary figure
figure('Name','Statistics','Renderer', 'painters', 'Position', [10 10 300 600])
clf
N       = model.n;
states  = cell2mat(simulatedM);
Nsample = size(states,2);
states  = zscore(states,0,2);
states  = reshape(states,N,N,3,Nsample);
fstates = fft(states,[],4);
pstates = reshape(mean(abs(fstates).^2,[1,2,3]),Nsample,1);
pstates = pstates(1:floor(Nsample/2));
pfreqs  = linspace(0,0.5/model.dt,floor(Nsample/2));
pstates = pstates./sum(pstates);
subplot(211);
plot(pfreqs(3:end).*60,pstates(3:end));
xlim([1,10]);
xlabel('Waves per minute');
ylabel('PSD (normalized)')
yticks([])
states  = cell2mat(simulatedM);
states  = reshape(states,N,N,3,Nsample);
q = reshape(states(:,:,1,:),N*N*Nsample,1);
a = reshape(states(:,:,2,:),N*N*Nsample,1);
r = reshape(states(:,:,3,:),N*N*Nsample,1);
subplot(212);
cla; hold all;
histogram(q,linspace(0,1,21),'Normalization','pdf','FaceColor',[0 0 1]);
histogram(a,linspace(0,1,21),'Normalization','pdf','FaceColor',[1 0 0]);
histogram(r,linspace(0,1,21),'Normalization','pdf','FaceColor',[0 1 0]);
xlim([0,1]);
xlabel('Occupancy fraction')
ylabel('Density (normalized)')
yticks([])
if savefigs, saveas(gcf,'Sampled_statistics_A.svg'); end





