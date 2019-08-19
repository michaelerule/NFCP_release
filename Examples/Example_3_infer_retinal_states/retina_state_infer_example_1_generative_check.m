%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up workspace
close all; clear all;
addpath('./','../','../../','../../NFCP'); 
NFCP_init

% Set this to the path of NFCP_datasets on your system
% (Contact mrule7404@gmail.com for a copy of NFCP_datasets)
datadir  = '/home/mer49/Dropbox (Cambridge University)/NFCP_datasets/';

% Set options

% Select dataset and binning resolution (10, 15 and 20 are available)
df       = 'P11_binned_100ms_20150219_leptin_retina1_control';
SRES     = 10;
subdir   = 'binned_spikes/spatiotemporal_binned/renewal_model/version_4_selected_datasets/';
fn       = sprintf('%d/%s_resolution_%d',SRES,df,SRES);
filepath = [datadir subdir fn '_preinitialized_using_renewal_model.mat'];
loadarchive 


% Heuristic model configuration
% Tuned for P11_binned_100ms_20150219_leptin_retina1_control
save_figures        = false;
skip_every          = 10; 
model.dt            = 0.1;%dt;  % time interval between observations
model.n             = 20;   % Simulation grid size
model.rAA           = 15;
model.sigma         = 0.1; % Standard-deviation for excitatory interaction kernel
model.thr           = 0.15;
model.update        = 'Laplace-subspace'; % Measurement update method
model.minrate       = 1e-6;
model.cutoff        = true; % Remove scales finer than interaction radius?
model.reg_state_var = 1e-4;  % Added uncertainty in state 
model.reg_count_var = 0;     % Added uncertainty in neuron count
model.reg_diag      = 1e-6;  % Added diagonal regularizaer
model.reg_inverse   = 1e-6;  % Diagonal regularizer for matrix inversions
model.ini_state_var = 1e-0;  % Initial variance in state estimate
model.ini_count_var = 0;     % Initial variance for neuron count
model.ini_reg_diag  = 1e-6;  % Initial variance added to diagonal
model.reg_precision = 1e-6;  % Regularization for precision matrix
model.alpha         = 1;
model.maxiter       = 25;
effectivepopsize    = 50000;
model.gain          = 100; % Gain for linear Cox-process observation model for spikes
model.bias          = 0;  % Bias for linear Cox-process observation model for spikes
model.alpha         = 1;  % Dispersion paramter, 1=Poisson
model.cutoff        = false;
model.verbosity     = 0;
model.safety        = 0;
model.cscale        = [1 3 1];
model.colors        = [0 1 0; 1 0 0; 0 0 1]; 
model.inhomogeneous_gain = 1;
model.names = strsplit('Q A R');
model.description = [ ...
%    Q  A  R rate
    -1  1  0 10    % spontaneous excitation 
     0 -1  1 2.5   % enter slow refractory loop
     1  0 -1 0.1   % slow refractory recovery
     ];

% Set system size rescaling based on effective population size
model = initializeModel(model);
model.ss_rescale    = 1/(model.dx*effectivepopsize);
fprintf('model.ss_rescale = %f\n',model.ss_rescale);

% Soft boundary around good regions
%U = goodchannels(:);%
U = ones(model.nn,1);
% U = gaussian2DblurOperator(model.n,0.1*model.n)*U;
ini = [U*0.8 U*0 U*0.1];

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
    'save_figure'     ,false
    'efraction'       ,0.04
    });


%{
fprintf('Animating...\n')
[llsum,infstate,margvar,infe] = stateInfer(ini,model,xydata,false,{...
    'doplot'       ,true 
    'upscale'      ,8
    'skipinf'      ,skip_every
    'showduration' ,1000  
    'showprog'     ,true  
    'showmaxy'     ,1
    'ratescale'    ,ratescale
    'peakactivity' ,false
    'points'       ,true
    'save_figure'  ,save_figures
    });
%}

