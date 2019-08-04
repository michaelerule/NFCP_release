%{
High-resolution sampling from Langevin model for illustration
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up workspace
close all; clear all; % Clean up namespace
rng('default'); 
rng('shuffle');      
addpath('./','../','../../','../../NFCP'); % Add paths
% Initialize library (adds more paths)
NFCP_init
% Start profiler
profile off; profile clear; profile on;  
% macros helpful for objective functions
functional_macros;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define model

% Number of "neurons" per unit area
% (simulation grid is defined as one unit area)
density = 50;

% Excitation model options
model.rQA     = 2e-1;  % Rate of spontaneous excitation (modeled as shot noise)
model.rAA     = 1.40./density;  % Excito-Excitatory interaction strength
model.thr     = 8e-3;  % Nonzero threshold for depolarization
model.sigma   = 0.1;   % Standard-deviation for excitatory interaction kernel

% Observation model options
model.gain    = 15; % Gain for linear Cox-process observation model for spikes
model.bias    = 0; % Bias for linear Cox-process observation model for spikes
model.alpha   = 1; % Dispersion paramter, 1=poisson

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

%{
% Slow refractory chain model
model.names = strsplit('Q A R1 R2');
model.description = [ ...
%    Q  A  R1  R2  rate
    -1  1   0   0   model.rQA % spontaneous excitation (zero because handled by rQA)
     0 -1   1   0   4e-1  % slow refractory loop
     0  0  -1   1   64e-4
     1  0   0  -1   64e-4 % slow refractory recovery
     ];
model.colors = [...
% Red  Green Blue (in range 0 to 1)
    0    1   0 
    1    0   0
    0    0   1
    0.3  0   0.7
    ];
%}

%{
% Slow-fast refractory model
model.names = strsplit('Q A R1 R2');
model.description = [ ...
%    Q  A  R1  R2  rate
    -1  1   0   0  model.rQA    % spontaneous excitation (zero because handled by rQA)
     0 -1   0   1  3e-2 % slow refractory loop
     1  0   0  -1  3e-4 % slow refractory recovery
     0 -1   1   0  3e-1 % fast refractory loop
     1  0  -1   0  3e-3 % fast refractory recovery
     ];
model.colors = [...
% Red  Green Blue (in range 0 to 1)
    0    1   0 
    1    0   0
    0    0   1
    0.3  0   0.7
    ];
%}

% Slow refractory model
model.names = strsplit('Q A R');
model.description = [ ...
%    Q  A  R1   rate
    -1  1   0   model.rQA     % spontaneous excitation (zero because handled by rQA)
     0 -1   1   4e-1  % slow refractory loop
     1  0  -1   32e-4 % slow refractory recovery
     ];
model.colors = [...
% Red  Green Blue (in range 0 to 1)
    0    1   0 
    1    0   0
    0    0   1
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select details of simulation and inference methods

% Space and time discretization and integration options
model.dt        = 1.0;   % time interval between observations
model.n         = 20;    % Simulation grid resolution
model.verbosity = 0;
model.safety    = 0;
mode.cutoff     = false; % if false, uses spatial grid size for population size
model.sigma     = 1.5/model.n;

% Plotting configuration
% Color scales for mapping Q/A/R fields for display
A_color_scale = 20.0;
Q_color_scale = 1.0 ;
R_color_scale = 1.0 ;
model.cscale = [Q_color_scale A_color_scale R_color_scale];
model = initializeModel(model);
[ini,xydata,rates,simulatedM] = stateSample(model,{...
    'doplot'      ,true 
    'Nsample'     ,1000 
    'Nburn'       ,250  
    'upscale'     ,8    
    'skipsim'     ,1
    'oversample'  ,10   
    'density'     ,density
    'save_figure' ,true
    });

