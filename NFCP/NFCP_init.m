%{
Initialization code for the NFCP library
%}
[basepath,name,ext] = fileparts(mfilename('fullpath'));
addpath([basepath '/util']);
addpath([basepath '/matrix']);
addpath([basepath '/output']);
addpath([basepath '/optimize']);
addpath([basepath '/spikes']);
addpath([basepath '/measurement']);

% These were removed 20190 08 18
% addpath([basepath '/EP']);
% addpath([basepath '/optimize/gpucb_gpml']);

functional_macros;
rng('shuffle');  
