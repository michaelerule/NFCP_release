%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gain adjustment implied by observation/prediction mismatch
K = model.nstates;      % number of species
N = model.n;            % side of NxN spatial grid
M = model.nn;           % number of spatial basis elements
T = size(xydata,2);
%{
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

% Go the other direction, to counts
Cest = bsxfun(@plus,bsxfun(@times,A,netgain),adjbias).*volume;
Cestok = Cest(ok,1:end);
Cok = counts(ok,1:end);
%}  
%{
% scatter(Cestok(:),Cok(:))
bins = prctile(Cestok(:),linspace(0,100,51))
for i=1:50
    a = bins(i);
    b = bins(i+1);
    inbin = (Cestok>=a)&(Cestok<b);
    x = Cestok(inbin);
    y = Cok(inbin);
    m(i) = mean(x)./model.ss_rescale;
    v(i) = var(y)./(model.ss_rescale.^2);
end
clf;
loglog(m,v);
factor = mean(v./m)
%}

%model.ss_rescale = factor;



states  = reshape(cell2mat(infstate),N,N,K,T);
total   = sum(states,3);
total   = reshape(total,N*N,T);
totalok = total(ok,T/3:T);
meantot = mean(total(:));

fprintf('System size should be rescaled by %f\n',meantot);

newss = model.ss_rescale/meantot
model.ss_rescale = newss;

newgain = model.gain*meantot
model.gain = newgain;

% Don't increase excitation
% on one hand, less noise means a smaller second-order correction
% to the excitation interaction, which we should compensate for
% on the other hand, a larger system size means smaller effective
% variance relative to the mean, which means that measurement updates
% have less effect as the precision of latent state estimation increses.
% This limits the ability of state measuremensts to counteract
% runaway excitation. 
%exc   = cell2mat(infe);
%maxRe = 0.5./(dt*prctile(exc(:),90))
%newexc = min(maxRe,model.rAA)
%mode.rAA = newexc;




