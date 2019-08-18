%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up workspace
close all; clear all; 
rng('default'); 
rng('shuffle');      
addpath('./','../','../../','../../NFCP'); 
NFCP_init
functional_macros;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATASET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n\n\n\n\n');
fprintf('______________________________________________________________________\n\n');
fprintf('Testing parameter estimates fit using a renewal model\n');
fprintf('______________________________________________________________________\n\n');

%{
directory = './20181103_renewal_model_processed/';

files = {
    'P10_18_Feb_15_leptin_ret1_ctl_SpkTs_dithered_binned_0.1s_preinitialized_using_renewal_model.mat',
    'P10 18 Feb 15 leptin ret2 ctl_SpkTs_dithered_binned_0.1s_preinitialized_using_renewal_model.mat',
    'P10_18_Feb_15_leptin_ret2_ctl_SpkTs_dithered_binned_0.1s_preinitialized_using_renewal_model.mat',
    'P10_6Feb15_ret1_leptin_ct_SpkTs_dithered_binned_0.1s_preinitialized_using_renewal_model.mat',
    'P10_6Feb15_ret2_leptin_ctl_SpkTs_dithered_binned_0.1s_preinitialized_using_renewal_model.mat',
    'P11_19_Feb_15_leptin_ret1_ctl_SpkTs_dithered_binned_0.1s_preinitialized_using_renewal_model.mat',
    'P11_19_Feb15_leptin_ret2_ctl_SpkTs_dithered_binned_0.1s_preinitialized_using_renewal_model.mat',
    'P11_25June12_day1_5_SpkTs_dithered_binned_0.1s_preinitialized_using_renewal_model.mat',
    'P4_18June12_Day2_4_SpkTs_dithered_binned_0.1s_preinitialized_using_renewal_model.mat',
    'P6_20June12_Day1_3_SpkTs_dithered_binned_0.1s_preinitialized_using_renewal_model.mat',
    'P6_25sept_spont1_SpkTs_dithered_binned_0.1s_preinitialized_using_renewal_model.mat',
    'P8_16Feb15_ret1_leptin_control_SpkTs_dithered_binned_0.1s_preinitialized_using_renewal_model.mat',
    'P8_16Feb15_ret2_leptin_control_SpkTs_dithered_binned_0.1s_preinitialized_using_renewal_model.mat'};
%}

% Want to test file
% 'P10_18_Feb_15_leptin_ret2_ctl_SpkTs_dithered_binned_0.1s_preinitialized_using_renewal_model.mat',
% output from revised model init script


N = 20;
directory = sprintf('./20181104/%d/',N);
files ={
	'P10_18_Feb_15_leptin_ret2_ctl_SpkTs_dithered_binned_0.1s_resolution_%d_preinitialized_using_renewal_model.mat',
	'P10_6Feb15_ret1_leptin_ct_SpkTs_dithered_binned_0.1s_resolution_%d_preinitialized_using_renewal_model.mat',
	'P11_25June12_day1_5_SpkTs_dithered_binned_0.1s_resolution_%d_preinitialized_using_renewal_model.mat',
	'P4_18June12_Day2_4_SpkTs_dithered_binned_0.1s_resolution_%d_preinitialized_using_renewal_model.mat',
	'P6_20June12_Day1_3_SpkTs_dithered_binned_0.1s_resolution_%d_preinitialized_using_renewal_model.mat',
	'P6_25sept_spont1_SpkTs_dithered_binned_0.1s_resolution_%d_preinitialized_using_renewal_model.mat'
	};


NFILES   = length(files);
fprintf('%d files prepared\n',NFILES);

FILENO   = 1;
fn       = sprintf(files{FILENO},N);

filepath = [directory fn];
fprintf('Using file %d at\n  %s\n  %s\n',FILENO,directory,fn);

loadarchive;

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

%{

average bias 7.563563
average gain 18721.082377

average bias 4.716316
average gain 116873.555569

average bias 75.461056
average gain 18699.768891
%}



% Initialize a basic model to call the binCounts function
model.dt = dt;  % time interval between observations
model.n  = N;   % Simulation grid size
model.interpolate = true; % interpolate when binning spikes:
NFRAMES  = size(xydata,2);


% Need to fix these!
ARRAY_DIMENSION_MM = 2.5
ARRAY_SIZE_MM2     = ARRAY_DIMENSION_MM^2

%{
% Convert from units of spikes/s∙mm² to spikes/bin_x²∙bin_t
bias = bias * ARRAY_SIZE_MM2 / (N*N) / dt;
% Convert from units of spikes/s∙mm² to spikes/ s∙array_x²
gain = gain * ARRAY_SIZE_MM2;

% CONVERT ROW MAJOR TO COLUMN MAJOR ORDER 
bias = reshape(bias,N,N)';
bias = bias(:);
inhomogeneous_gain = reshape(inhomogeneous_gain,N,N)';
inhomogeneous_gain = inhomogeneous_gain(:);
inhomogeneous_gain(bias==0) = 0;
%}

gain = gain./0.7;

fprintf('gain=%f\n',gain);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Observation model options
model.gain = gain;%/0.8;
model.bias = bias;
model.inhomogeneous_gain = inhomogeneous_gain;
model.alpha  = 1;%./alpha;
model.volume = dt/N.^2;
model.ss_rescale = 1;
model.rAA = rE;%*3.5;   % Excito-Excitatory interaction strength

initialize_five_state_model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fix the gain model

figure('name','bias/gain sanity check')
binned = {};
for i=1:size(xydata,2)
    binned{i} = binCounts(model,xydata{1,i}); 
end
binned  = cell2mat(binned);
counts  = binned ./ model.alpha;
% Counts are in spikes/bin_x²∙bin_t
% Convert to unitless
fraction = bsxfun(@rdivide,counts,(model.volume.*model.gain.*model.inhomogeneous_gain));
fraction = bsxfun(@minus,fraction,model.adjusted_bias);
%fraction = fraction(model.bias>0,1:end);


% Sanity check that gain/bias arrays are not transposed
subplot(221)
imagesc(reshape(mean(counts,2),N,N));
title('mean count');
subplot(222);
imagesc(reshape(bias,N,N));
title('bias');
subplot(223);
imagesc(reshape(inhomogeneous_gain,N,N));
title('gain');

fprintf('average bias %f\n',mean(model.bias))
fprintf('average gain %f\n',model.gain)

animatestate


