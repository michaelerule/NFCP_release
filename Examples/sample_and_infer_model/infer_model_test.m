% Set up workspace
close all; clear all; rng('shuffle');      % Clean up namespace
addpath('./','../','../../','../../NFCP'); % Add library path (relative path)
NFCP_init

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define model

% Excitation model options
model.rQA     = 2e-1;  % Rate of spontaneous excitation (modeled as shot noise)
density       = 50;    % Number of "neurons" per unit area
model.rAA     = 1.40./density;  % Excito-Excitatory interaction strength
model.thr     = 8e-3;  % Nonzero threshold for depolarization
model.sigma   = 0.1;   % Standard-deviation for excitatory interaction kernel

% Observation model options
model.gain    = 10; % Gain for linear Cox-process observation model for spikes
model.bias    = 0;  % Bias for linear Cox-process observation model for spikes
model.alpha   = 1;  % Dispersion paramter, 1=poisson

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define linear (spontaneous) transitions
% Format is:
%   Matrix
%   nrows = number of transitions
%   ncols = number of states plus one (extra column for the rates)
%   row,column = change in state for this reaction
%   row,colum = -1 --> neuron leaves this state
%   row,colum =  1 --> neuron enters this state
%   Only one neuron should change states for each reaction
%   Last column: reaction rates, in particles / time unit

% Slow refractory model
model.names = strsplit('Q A R');
%model.description = [ ...
%  Q  A  R1   rate
% -1  1   0   model.rQA     % spontaneous excitation (zero because handled by rQA)
%  0 -1   1   4e-1  % slow refractory loop
%  1  0  -1   32e-4 % slow refractory recovery
% ];
     
model.rAA = 1.4./density;  % Excito-Excitatory interaction strength
model.description = [ ...
% Q  A R1   rate
 -1  1  0   model.rQA     % spontaneous excitation (zero because handled by rQA)
  0 -1  1   5e-1 % slow refractory loop
  1  0 -1   1e-2 % slow refractory recovery
  ];

model.colors = [...
% R G B
  0 1 0 
  1 0 0
  0 0 1
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select details of simulation and inference methods

% Space and time discretization and integration options
model.dt        = 1.0;   % time interval between observations
model.n         = 8;     % Simulation grid resolution, i.e. n^2 basis functionss
model.cutoff    = false; % Use interaction scales to set minimum spatial scale?
model.safety    = 0
model.verbosity = 0
model.bias      = 1e-9;
model.link      = 'linear'
% Plotting configuration
% Color scales for mapping Q/A/R fields for display
model.cscale = [1 20 1]; % Q A R
[ini,xydata,rates,simulatedM] = stateSample(model,{...
    'doplot'      ,true 
    'Nsample'     ,1000 
    'Nburn'       ,250  
    'upscale'     ,8    
    'skipsim'     ,100
    'oversample'  ,20   
    'density'     ,density
    'save_figure' ,false
    });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjustments for state-inference
% Turn off the spontaneous exictation (infer as extrinsic noise)
profile off; profile clear; profile on;

model.update         = 'Laplace'; % Measurement update method
model.rQA            = 0.0;
model.maxiter        = 500
model.linearRates(1) = 0.0;
model.minrate        = 1e-4;
model.cutoff         = false
model.reg_state_var  = 1e-8;
model.reg_count_var  = 0;
model.reg_diag       = 1e-8;
model.reg_inverse    = 1e-12;
model = initializeModel(model);
[llsum,infstate,margvar,infe] = stateInfer(ini,model,xydata,simulatedM,{...
    'doplot'       ,true 
    'upscale'      ,8
    'skipinf'      ,100
    'showduration' ,1000  
    'showmaxy'     ,density+5   
    'ratescale'    ,50
    'peakactivity' ,false
    'points'       ,false
    });
profile viewer;

%coordinate_sweep;
