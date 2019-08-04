%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up workspace
close all; clear all;
addpath('./','../','../../','../../NFCP'); 
NFCP_init

directory = './';
%fn = 'P10_binned_100ms_20150218_leptin_retina2_control_resolution_10_preinitialized_using_renewal_model';
%fn = 'P10_binned_100ms_20150218_leptin_retina2_control_resolution_20_preinitialized_using_renewal_model';
fn = 'P11_binned_100ms_20120514_spontaneous_day2_3_resolution_10_preinitialized_using_renewal_model.mat';
%fn = 'P11_binned_100ms_20120625_day1_5_resolution_20_preinitialized_using_renewal_model';
filepath = [directory fn];
loadarchive

model.dt      = dt;  % time interval between observations
model.n       = N;   % Simulation grid size
model.rAA     = 20;  % Excito-Excitatory interaction strength
model.sigma   = 0.15; % Standard-deviation for excitatory interaction kernel
model.gain    = 30000;% Gain for linear Cox-process observation model for spikes
model.bias    = 0;   % Bias for linear Cox-process observation model for spikes
model.alpha   = 1;   % Dispersion paramter, 1=poisson
model.interpolate = true;
model.inhomogeneous_gain = inhomogeneous_gain;
model.bias = bias/(N*N);
% Slow refractory model
model.names = strsplit('Q A R1 R2');
model.description = [ ...
%    Q  A R1 R2 rate
    -1  1  0  0 0    % spontaneous excitation 
     0 -1  1  0 5e0   % enter slow refractory loop
     0  0 -1  1 8e-1  % slow refractory recovery
     1  0  0 -1 8e-1  % slow refractory recovery
     ];
model.cscale        = [1.1 5 1.5 1.5];
model.colors        = [0 1 0; 1 0 0; 0 0 1; 0 0 1]; 
model.thr           = 1e-2;
model.update        = 'Laplace-subspace'; % Measurement update method
model.minrate       = 1e-6;
model.cutoff        = false; % Remove scales finer than interaction radius?
model.reg_state_var = 1e-6;  % Added uncertainty in state 
model.reg_count_var = 1e-9;  % Added uncertainty in neuron count
model.reg_diag      = 1e-9;  % Added diagonal regularizaer
model.reg_inverse   = 1e-9;  % Diagonal regularizer for matrix inversions
model.ini_state_var = 1e02;  % Initial variance in state estimate
model.ini_count_var = 1e-9;  % Initial variance for neuron count
model.ini_reg_diag  = 1e-9;  % Initial variance added to diagonal
model.reg_precision = 1e-9;  % Regularization for precision matrix

model.alpha = 1;
model.ss_rescale = 30;

model = initializeModel(model);


% Re-compute gains/biases
%{
+--------------------------+-------------------------------------------+
| Variable                 | Units                                     | 
+==========================+===========================================+
| model.bias               | spikes       /   bin_x²∙bin_t             |
+--------------------------+-------------------------------------------+
| model.gain               | spikes       /   s∙array_x²               |
+--------------------------+-------------------------------------------+
%}
counts = {};
for i=1:size(xydata,2),
    % Create a histogram of 2D xypoints (the observation counts)
    % These are premultiplied by alpha
    counts{i} = binCounts(model,xydata{i});  
end
counts = cell2mat(counts);
for i=1:N*N,
    b2(i) = prctile(counts(i,1:end),50);
end
for i=1:N*N,
    g2(i) = prctile(counts(i,1:end),99.9);
end
goodchannels = inhomogeneous_gain>0;
g2 = g2(:);
model.gain = mean(g2(goodchannels)./inhomogeneous_gain(goodchannels))/model.volume;
model.bias = b2;

fprintf('Old gain was %0.4f\n',gain);
fprintf('New gain is  %0.4f\n',model.gain);
fprintf('Old bias was %0.4f\n',mean(bias(goodchannels)));
fprintf('New bias is  %0.4f\n',mean(model.bias(goodchannels)));

% Soft boundary around good regions
%U = goodchannels(:);%
U = ones(model.nn,1);
% U = gaussian2DblurOperator(model.n,0.1*model.n)*U;
ini = [U*0.2 U*0 U*0.4 U*0.4];


fprintf('Animating...\n')
[llsum,infstate,margvar,infe] = stateInfer(ini,model,xydata,false,{...
    'doplot'       ,true 
    'upscale'      ,8
    'skipinf'      ,1
    'showduration' ,1000  
    'showmaxy'     ,1
    'ratescale'    ,5
    'peakactivity' ,false
    'points'       ,true
    'save_figure'  ,false
    });



