%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up workspace
close all; clear all; rng('shuffle');      
addpath('./','../','../../','../../NFCP'); % Add paths
NFCP_init

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define model

% Number of "neurons" per unit area
% (simulation grid is defined as one unit area)
effectivepopsize = 50;

% Excitation model options
model.rAA     = 1.40./effectivepopsize;  % Excito-Excitatory interaction strength
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
    'effectivepopsize'     ,effectivepopsize
    'save_figure' ,false
    });

%{
Returned variables are: 

ini : `matrix`, 
    Packed initial conditions for mean concentrations. Species are
    packed in the order defined by the model, where the second
    species reflects the active "A" state. 

xydata : `cell`, 1×Ntimes 
    Cell array of spiking events. For each timepoint, we have a 
    Nspikes×2 matrix, where column 1 is the x location and column 2 
    is the y location of each spike, in the [0,1]² box.

rates: `cell`, 1×Ntimes 
    Basis-projected spiking events, with rates normalized by
    spatiotemporal volume. 

simulatedM: `cell`, 1×Ntimes 
    Each cell array contains a MN² by 1 matrix of packed states, where 
    M is the number of states (species), and N is the size of the N×N 
    spatial basis, with N² basis elements total. Species are packed in 
    order, with the spatial values packed in Matlab's default 
    (column major) order. 
%}


%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Infer states while animating
%}


% Adjustments for state-inference
% Turn off the spontaneous exictation (infer as extrinsic noise)
model.rQA               = 0.0;
model.linearRates(1)    = model.rQA;
% Turn off finite threshold on means 
model.thr               = 0.0;
model.update            = 'Laplace-subspace'; % Measurement update method
model.cutoff            = false
infermodel.dolikelihood = false;
model = initializeModel(model);

fprintf('Animating...')
tic()
profile off; profile clear; profile on;  
[llsum,infstate,margvar,infe] = stateInfer(ini,model,xydata,simulatedM,{...
    'doplot'       ,true 
    'upscale'      ,8
    'skipinf'      ,100
    'showduration' ,1000  
    'showmaxy'     ,effectivepopsize+5   
    'ratescale'    ,25
    'peakactivity' ,false
    'points'       ,false
    'save_figure'  ,false
    });
profile off; profile viewer;
toc()


