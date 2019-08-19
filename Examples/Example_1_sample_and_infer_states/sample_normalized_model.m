%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up workspace
close all; clear all; rng('shuffle');      
addpath('./','../','../../','../../NFCP'); % Add paths
NFCP_init

%{
Example 1c: 

Sample from a normalized model and infer states
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define model

% Number of "neurons" per unit area (sim. grid is defined as one unit area)
effectivepopsize = 50;

% Excitation model options
model.rAA     = 1.40;  % Excito-Excitatory interaction strength
model.thr     = 8e-3;  % Nonzero threshold for depolarization
model.sigma   = 0.1;   % Standard-deviation for excitatory interaction kernel

% Observation model options
model.gain    = 15*effectivepopsize; % Gain for linear Cox-process observation model for spikes
model.bias    = 0;  % Bias for linear Cox-process observation model for spikes
model.alpha   = 1;  % Dispersion paramter, 1=Poisson

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
model.dt        = 1.0; % time interval between observations
model.n         = 12;  % Simulation grid resolution
model.verbosity = 0;
model.safety    = 0;
model.cscale    = [1 20 1]; % Color scales for mapping Q/A/R for display

model.cutoff    = false;
model.ss_rescale= 1./effectivepopsize;

model = initializeModel(model);
Nsample = 1000;
fprintf(1,'Sampling from model...\n');
[ini,xydata,rates,simulatedM] = stateSample(model,{...
    'doplot'      ,true 
    'Nsample'     ,Nsample 
    'Nburn'       ,250  
    'upscale'     ,8    
    'skipsim'     ,20
    'oversample'  ,10
    'effectivepopsize'     ,1
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
Infer states while animating, using a normalized model.
%}

infermodel = initializeModel(model);

% Adjustments for state-inference
infermodel.rQA            = 0.0;
infermodel.linearRates(1) = 0.0;
infermodel.thr            = 0.0;
infermodel.update         = 'Laplace-subspace'; 
infermodel.cutoff         = false;
infermodel.dolikelihood   = false;

% The state inference should converge from uncertain initial conditions, but it
% takes O(10000) timesteps to do so. For a quick illustration, we start the
% inference around the correct initial conditions.
infermodel.ini_state_var  = 1e-2;
infermodel.minrate        = 1e-4./effectivepopsize;

% By default the code will try to scale fluctuations by 1/N²
% In the normalized model, we want fluctuations to be scaled by 1/effectivepopsize.
% This factor overrides the default scaling to ensure this:
%if infermodel.cutoff,
%    infermodel.ss_rescale = 1./(effectivepopsize*2*pi*model.sigma^2.);
%else
%    infermodel.ss_rescale = model.n.^2./effectivepopsize;
%end

infermodel = initializeModel(infermodel);

fprintf('Animating...')
tic()
profile off; profile clear; profile on;  
[llsum,infstate,margvar,infe] = stateInfer(ini,infermodel,xydata,simulatedM,{...
    'doplot'       ,true 
    'upscale'      ,8
    'skipinf'      ,20
    'showduration' ,1000
    'showmaxy'     ,1.05   
    'ratescale'    ,25
    'peakactivity' ,false
    'points'       ,false
    'save_figure'  ,false
    'showprog'     ,true
    });
profile off; profile viewer;
toc()


