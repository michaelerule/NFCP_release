%{

There was a bug in the units for the observation model.
This script was used to correct that. 

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up workspace
close all; clear all; 
rng('default'); 
rng('shuffle');      
addpath('./','../','../../','../../NFCP'); 
NFCP_init
functional_macros;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILD MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

    +--------------------------+----------------+---------------------+
    | Variables                | Units...       | ...per units        | 
    +==========================+================+=====================+
    | model.dt                 | seconds      / | bin_t               |
    +--------------------------+----------------+---------------------+
    | model.dx                 | array_x²     / | bin_x²              |
    +--------------------------+----------------+---------------------+
    | model.n                  | bin_x        / | array_x             |
    +--------------------------+----------------+---------------------+
    | model.nn                 | bin_x²       / | array_x²            |
    +--------------------------+----------------+---------------------+
    | model.volume             | s∙array_x²   / | bin_x²∙bin_t        |
    +--------------------------+----------------+---------------------+
    | model.bias               | spikes       / | bin_x²∙bin_t        |
    +--------------------------+----------------+---------------------+
    | model.gain               | spikes       / | s∙array_x²          |
    +--------------------------+----------------+---------------------+
    | model.alpha              | spike mean   / | spike standard dev. |
    +--------------------------+----------------+---------------------+
    | model.inhomogeneous_gain | unitless                             |
    +--------------------------+----------------+---------------------+
    | model.adjusted_bias      | unitless                             |
    +--------------------------+----------------+---------------------+
    | model.adjusted_gain      | alpha∙spikes / | bin_x²∙bin_t        |
    +--------------------------+----------------+---------------------+
    | rate (inside obj. fun.)  | spikes       / | bin_x²∙bin_t        |
    +--------------------------+----------------+---------------------+

%}    



directory1 = './20181103_renewal_model_processed/';
directory2 = './'%20181104/10/';
fn1 = 'P10_18_Feb_15_leptin_ret2_ctl_SpkTs_dithered_binned_0.1s_preinitialized_using_renewal_model.mat'
fn2 = 'P10_18_Feb_15_leptin_ret2_ctl_SpkTs_dithered_binned_0.1s_resolution_10_preinitialized_using_renewal_model.mat',
	
fn = fn1
filepath = [directory1 fn1];
loadarchive;
oldgain = gain
oldbias = bias
% Initialize a basic model to call the binCounts function
model.dt = dt;  % time interval between observations
model.n  = N;   % Simulation grid size
model.interpolate = true; % interpolate when binning spikes:
NFRAMES  = size(xydata,2);
% Need to fix these!
ARRAY_DIMENSION_MM = 2.5
ARRAY_SIZE_MM2     = ARRAY_DIMENSION_MM^2
% Convert from units of spikes/s∙mm² to spikes/bin_x²∙bin_t
bias = bias .* (ARRAY_SIZE_MM2 / (N*N) / dt);
% Convert from units of spikes/s∙mm² to spikes/ s∙array_x²
gain = gain * ARRAY_SIZE_MM2;
% CONVERT ROW MAJOR TO COLUMN MAJOR ORDER 
bias = reshape(bias,N,N)';
bias = bias(:);
inhomogeneous_gain = reshape(inhomogeneous_gain,N,N)';
inhomogeneous_gain = inhomogeneous_gain(:);
inhomogeneous_gain(bias==0) = 0;
model.gain = gain;
model.bias = bias;
model.inhomogeneous_gain = inhomogeneous_gain;
model1 = model;

fn = fn2
filepath = [directory2 fn2];
loadarchive;
% Initialize a basic model to call the binCounts function
model.dt = dt;  % time interval between observations
model.n  = N;   % Simulation grid size
model.interpolate = true; % interpolate when binning spikes:
model.gain = gain;
model.bias = bias;
model.inhomogeneous_gain = inhomogeneous_gain;

oldgain
model1.gain
model.gain
% Units are completely off
mean(oldbias)
mean(model1.bias)
mean(model.bias)


