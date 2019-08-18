%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up workspace
close all; clear all; 
rng('default'); 
rng('shuffle');      
addpath('./','../','../../','../../NFCP'); 
NFCP_init

directory = './'
fn = 'P10_binned_100ms_20150218_leptin_retina2_control_resolution_10_preinitialized_using_renewal_model'
filepath = [directory fn];
loadarchive

%retina_state_infer;


