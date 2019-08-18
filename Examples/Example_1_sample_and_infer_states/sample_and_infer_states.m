%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up workspace
close all; clear all; rng('shuffle');      
addpath('./','../','../../','../../NFCP'); % Add paths
NFCP_init

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define model

% Number of "neurons" per unit area
% (simulation grid is defined as one unit area)
density = 50;

% Excitation model options
model.rAA     = 1.40./density;  % Excito-Excitatory interaction strength
model.thr     = 8e-3;  % Nonzero threshold for depolarization
model.sigma   = 0.1;   % Standard-deviation for excitatory interaction kernel

% Observation model options
model.gain    = 15; % Gain for linear Cox-process observation model for spikes
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

%{
% Slow refractory chain model
model.names = strsplit('Q A R1 R2');
model.description = [ ...
%    Q  A R1 R2  rate
    -1  1  0  0  2e-1  % spontaneous excitation
     0 -1  1  0  4e-1  % slow refractory loop
     0  0 -1  1  64e-4
     1  0  0 -1  64e-4 % slow refractory recovery
     ];
% RGB tuples for all species
model.colors = [0 1 0; 1 0 0; 0 0 1; 0.3 0 0.7];
%}

%{
% Slow-fast refractory model
model.names = strsplit('Q A R1 R2');
model.description = [ ...
%    Q  A R1 R2  rate
    -1  1  0  0  2e-1 % spontaneous excitation
     0 -1  0  1  3e-2 % slow refractory loop
     1  0  0 -1  3e-4 % slow refractory recovery
     0 -1  1  0  3e-1 % fast refractory loop
     1  0 -1  0  3e-3 % fast refractory recovery
     ];
% RGB tuples for all species
model.colors = [0 1 0; 1 0 0; 0 0 1; 0.3 0 0.7];
%}

% Slow refractory model
model.names = strsplit('Q A R');
model.description = [ ...
%    Q  A R1  rate
    -1  1  0  2e-1  % spontaneous excitation 
     0 -1  1  4e-1  % slow refractory loop
     1  0 -1  32e-4 % slow refractory recovery
     ];

% RGB tuples for all species
model.colors = [0 1 0; 1 0 0; 0 0 1]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select details of simulation and inference methods

% Space and time discretization and integration options
model.dt        = 1.0;           % time interval between observations
model.n         = 10;             % Simulation grid resolution
model.verbosity = 0;
model.safety    = 0;
model.cscale    = [1 20 1]; % Color scales for mapping Q/A/R for display
model = initializeModel(model);
[ini,xydata,rates,simulatedM] = stateSample(model,{...
    'doplot'      ,true 
    'Nsample'     ,1000 
    'Nburn'       ,250  
    'upscale'     ,8    
    'skipsim'     ,20
    'oversample'  ,10   
    'density'     ,density
    'save_figure' ,false
    });

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Infer states while animating
%}


% Adjustments for state-inference
% Turn off the spontaneous exictation (infer as extrinsic noise)
model.rQA            = 0.0;
model.linearRates(1) = 0.0;
% Turn off finite threshold on means 
% (this would be properly handled as truncated normal)
model.thr            = 0.0;
model.update         = 'Laplace'; % Measurement update method
model.cutoff         = false
model = initializeModel(model);

fprintf('Animating...')
tic()
profile off; profile clear; profile on;  
[llsum,infstate,margvar,infe] = stateInfer(ini,model,xydata,simulatedM,{...
    'doplot'       ,true 
    'upscale'      ,8
    'skipinf'      ,10
    'showduration' ,1000  
    'showmaxy'     ,density+5   
    'ratescale'    ,25
    'peakactivity' ,false
    'points'       ,false
    'save_figure'  ,false
    });
profile off; profile viewer;
toc()


