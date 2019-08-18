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
for i=1:T
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
states = reshape(cell2mat(infstate),N,N,K,T);
A      = reshape(states(1:N,1:N,2,1:T),N*N,T);
ok         = model.goodchannels;
Aok        = A(ok,1:end);
Aestok     = Aest(ok,1:end);
model.gain = model.gain * mean(Aestok(:)) / mean(Aok(:));


