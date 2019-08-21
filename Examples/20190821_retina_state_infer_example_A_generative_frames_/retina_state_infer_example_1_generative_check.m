%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up workspace
close all; clear all;
addpath('./','../','../../','../NFCP'); 
NFCP_init

% Set this to the path of NFCP_datasets on your system
% (Contact mrule7404@gmail.com for a copy of NFCP_datasets)
datadir  = '/home/mer49/Dropbox (Cambridge University)/NFCP_datasets/';

% Set options

% Select dataset and binning resolution (10, 15 and 20 are available)
df       = 'P11_binned_100ms_20150219_leptin_retina1_control';
SRES     = 20;
subdir   = 'binned_spikes/spatiotemporal_binned/renewal_model/version_4_selected_datasets/';
fn       = sprintf('%d/%s_resolution_%d',SRES,df,SRES);
filepath = [datadir subdir fn '_preinitialized_using_renewal_model.mat'];
loadarchive 

% Heuristic model configuration
% Tuned for P11_binned_100ms_20150219_leptin_retina1_control
save_figures        = false;
skip_every          = 10; 
model.dt            = 0.1;   % time interval between observations
model.n             = SRES;    % Simulation grid size
model.rAA           = 15;
model.sigma         = 0.05;   % Standard-deviation for excitatory interaction kernel
model.thr           = 0.15;
model.update        = 'Laplace-subspace'; % Measurement update method
model.minrate       = 1e-6;
model.cutoff        = true;  % Remove scales finer than interaction radius?
model.reg_state_var = 1e-4;  % Added uncertainty in state 
model.reg_count_var = 0;     % Added uncertainty in neuron count
model.reg_diag      = 1e-6;  % Added diagonal regularizaer
model.reg_inverse   = 1e-6;  % Diagonal regularizer for matrix inversions
model.ini_state_var = 1e-0;  % Initial variance in state estimate
model.ini_count_var = 0;     % Initial variance for neuron count
model.ini_reg_diag  = 1e-6;  % Initial variance added to diagonal
model.reg_precision = 1e-6;  % Regularization for precision matrix
model.maxiter       = 25;
effectivepopsize    = 100;
model.gain          = 0;     % Gain for linear Cox-process observation model for spikes
model.bias          = 0;     % Bias for linear Cox-process observation model for spikes
model.alpha         = 1;     % Dispersion paramter, 1=Poisson
model.cutoff        = false;
model.verbosity     = 0;
model.safety        = 0;
model.cscale        = [1 3 1];
model.colors        = [0 0 1; 1 0 0; 0 1 0];
model.inhomogeneous_gain = 1;
model.names = strsplit('Q A R');
model.description = [ ...
%    Q  A  R rate
    -1  1  0 10    % spontaneous excitation 
     0 -1  1 1.5   % enter slow refractory loop
     1  0 -1 0.1   % slow refractory recovery
     ];

% Set system size rescaling based on effective population size
model = initializeModel(model);
model.ss_rescale = 1/effectivepopsize;
fprintf('model.ss_rescale = %f\n',model.ss_rescale);

% Soft boundary around active regions
% (messy code gets reflected boundary conditions for a spatial blur)counts = {};
for i=1:size(xydata,2), 
    counts{i} = binCounts(model,xydata{i}); 
end
counts = cell2mat(counts);
for i=1:N*N, 
    m2(i) = mean(counts(i,1:end)); 
end
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

ratescale = 5;
model.names{2} = sprintf('A×%d',ratescale);

model = initializeModel(model);
Nsample = 1000;

fprintf(1,'Sampling from model...\n');
[ini,xydata,rates,simulatedM] = stateSample(model,{...
    'doplot'          ,true 
    'Nsample'         ,Nsample 
    'Nburn'           ,250  
    'upscale'         ,8    
    'skipsim'         ,skip_every
    'oversample'      ,10   
    'effectivepopsize',1
    'save_figure'     ,true
    'efraction'       ,0.04
    'ini'             ,ini
    'rescale'         ,true
    });



