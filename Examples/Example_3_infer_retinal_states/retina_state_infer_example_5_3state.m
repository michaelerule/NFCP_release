%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up workspace
close all; clear all;
addpath('./','../','../../','../../NFCP'); 
NFCP_init

directory = './';
fn = 'P11_binned_100ms_20150219_leptin_retina1_control_resolution_20_preinitialized_using_renewal_model.mat';
filepath = [directory fn];
loadarchive

model.dt      = dt;  % time interval between observations
model.n       = N;   % Simulation grid size
model.rAA     = 20;   % Excito-Excitatory interaction strength
model.sigma   = 0.1; % Standard-deviation for excitatory interaction kernel
model.inhomogeneous_gain = inhomogeneous_gain;
% Slow refractory model
model.names = strsplit('Q A R');
model.description = [ ...
%    Q  A  R rate
    -1  1  0 0     % spontaneous excitation 
     0 -1  1 20    % enter slow refractory loop
     1  0 -1 0.07  % slow refractory recovery
     ];
model.cscale        = [1 3 1];
model.colors        = [0 1 0; 1 0 0; 0 0 1]; 
model.thr           = 0.1;
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
model = initializeModel(model);

effectivepopsize    = 1;
model.ss_rescale    = 1/(model.dx*effectivepopsize);
fprintf('model.ss_rescale = %f\n',model.ss_rescale);
model.maxiter = 25;

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
    g2(i) = prctile(counts(i,1:end),99.99);
end
goodchannels = inhomogeneous_gain>0;
g2 = g2(:);
model.gain = mean(g2(goodchannels)./inhomogeneous_gain(goodchannels))/model.volume;
model.bias = b2;

%model.bias = bias/(N*N);
%model.gain = model.gain/model.alpha;

model.inhomogeneous_gain = gaussian2DblurOperator(model.n,model.sigma.*model.n)*model.inhomogeneous_gain;
model.bias = gaussian2DblurOperator(model.n,model.sigma.*model.n)*model.bias';

fprintf('Old gain was %0.4f\n',gain);
fprintf('New gain is  %0.4f\n',model.gain);
fprintf('Old bias was %0.4f\n',mean(bias(goodchannels)));
fprintf('New bias is  %0.4f\n',mean(model.bias(goodchannels)));

% Soft boundary around good regions
%U = goodchannels(:);%
U = ones(model.nn,1);
% U = gaussian2DblurOperator(model.n,0.1*model.n)*U;
ini = [U*0.8 U*0 U*0.1];

ratescale = 5;
model.names{2} = sprintf('A×%d',ratescale);
fprintf('Animating...\n')
[llsum,infstate,margvar,infe] = stateInfer(ini,model,xydata,false,{...
    'doplot'       ,true 
    'upscale'      ,8
    'skipinf'      ,1
    'showduration' ,1000  
    'showprog'     ,true  
    'showmaxy'     ,1
    'ratescale'    ,ratescale
    'peakactivity' ,false
    'points'       ,true
    'save_figure'  ,true
    });

%{
[llsum,infstate,margvar,infe] = stateInfer(ini,model,xydata(1:500),false,{...
    'doplot'       ,false 
    'showprog' ,true  
    });


%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gain adjustment implied by observation/prediction mismatch
K = model.nstates;      % number of species
N = model.n;            % side of NxN spatial grid
M = model.nn;           % number of spatial basis elements
T = size(xydata,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fix the gain model
binned = {};
for i=1:500%T
    binned{i} = binCounts(model,xydata{1,i}); 
end
binned  = cell2mat(binned);
counts  = binned ./ model.alpha;
volume  = model.volume;
netgain = model.gain * model.inhomogeneous_gain(:);
adjbias = model.bias(:).^model.gamma;
Aest    = bsxfun(@rdivide,bsxfun(@minus,counts/volume,adjbias),netgain);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gain adjustment implied by observation/prediction mismatch
% Get inferred states as a n x n x nstated x ntimes array
states = reshape(cell2mat(infstate),N,N,K,500);%);T);
%A      = reshape(states(1:N,1:N,2,1:T),N*N,500);%T);
A      = reshape(states(1:N,1:N,2,1:500),N*N,500);%T);
ok         = model.goodchannels;
Aok        = A(ok,1:end);
Aestok     = Aest(ok,1:end);
model.gain = model.gain * mean(Aestok(:)) / mean(Aok(:));


%}
