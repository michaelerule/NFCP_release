% Parameter corrections

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Don't use all time points
Nskip     = 100;

K = model.nstates;      % number of species
N = model.n;            % side of NxN spatial grid
M = model.nn;           % number of spatial basis elements
T = size(xydata,2);

% Rate of change of intensities for the Q, A, R fields
x    = cell2mat(infstate);
dxdt = diff(x,1,2)./model.dt;

% Factors affecting rate of change
% The state intensities Q, A, R determine their rates of change
% The spontaneous Q->A state is disabled, and this transition
% instead dpends on the effective excitation (`excitation`)
% So we overwrite the Q concentration with effective excitation
y = x(1:end,1:T-1); % copy, exclude last timepoint
% N²xT excitatory inputs (inferred during filtering)
Q = y(1:M,1:end);
A = y(M+1:M*2,1:end);
y(1:M,1:end) = cell2mat(infe(1,1:T-1));



% Instead compute the active state directly from the counts?
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
fraction = max(1e-6,fraction);
fraction = fraction(1:end,1:T-1);
%fraction = fraction(model.bias>0,1:end);
%overwrite A state?
%A = fraction*0.1+A*0.9;
y(M+1:M*2,1:end) = A;

% We should recompute <QA>
%model.K2D = gaussian2DblurOperator(model.n,model.sigma*model.n);
%model.K1D = gaussian1DblurOperator(model.n,model.sigma*model.n);
%model.blur = @(x) reshape(model.K1D*reshape(x,model.n,model.n)*model.K1D',1,model.nn);
y(1:M,1:end) = Q.*(model.K2D*A);

% Shape data as Nbases x Nspecies x Time-1
% X is the change in intensities that we will try to predict
% Y contains the variables upon which dX/dT depends
x_spatial = reshape(dxdt,[M K T-1]); % Nspexies x Ntimes*Nspatial matrix
y_spatial = reshape(y   ,[M K T-1]); % Nspexies x Ntimes*Nspatial matrix

x_spatial = x_spatial(1:M,1:K,1:Nskip:T-1);
y_spatial = y_spatial(1:M,1:K,1:Nskip:T-1);
Tkept     = size(x_spatial,3);

% Collapse space/time into single dimension
% First, reorder data as Nspecies x Nbases x Time
% Then collapse the Nbases and Time dimensions
% Forming a single Nspecies x Ndata array
% TODO: consider removing or masking "bad" channels
x = reshape(permute(x_spatial,[2 1 3]),[K M*Tkept]);
y = reshape(permute(y_spatial,[2 1 3]),[K M*Tkept]);

% Set up optimization problem

% Use true values as initial conditions
% We keep all refractory transitions at the same rate, so there will
% be fewer rate parameters than species
Nrefract = K-2;
% Use true values as initial conditions
Re = model.rAA;
Ra = model.linearRates(2);
Rr = model.linearRates(3);
H0 = [Re,Ra,Rr];
Nparams  = numel(H0);
parnames = {'Excitation rate','A→R rate','R→R and R→Q rate'};

% Run a test example
% Predicted transition rates are
% Parameters * Regressors (concentrations, Y)
% But we collapsed all refractor rates into one parameter!
% We need to expand it
E = [1 0 0 0 0
     0 1 0 0 0
     0 0 1 1 1]';

% Express problem using lambda calculus
dXhat     = @(RY) RY([5 1 2 3 4])-RY; % reactions
loss      = @(R,dX,Y) norm(dXhat((E*R').*Y)-dX).^2; % squared loss
objective = @(H,dX,Y) loss(exp(H),dX,Y); % objective

options = optimset('MaxIter',1e12,'MaxFunEvals',1e12,'TolFun',1e-4);

% Solve many small local optimization problems
% Subsample to avoid doing more work than necessary
L       = size(x,2);
result  = {};
parfor_progress(L);
for i=1:L
    Y   = y(1:end,i);
    dX  = x(1:end,i);
    fun = @(H) objective(H,dX,Y);
    result{i} = exp(fminsearch(fun,log(H0),options))';
    parfor_progress;
end
result = cell2mat(result);
parfor_progress(0);

result_reshaped = reshape(result,Nparams,N,N,Tkept);

ok = model.goodchannels;

% Get the mean at each (valid/good) spatial location
rmeans = nanmean(result_reshaped,3);
for i=1:Nparams,
    parameter_ests   = rmeans(i,ok);
    parameter_mean   = mean(parameter_ests(:));
    parameter_median = median(parameter_ests(:));
    parameter_var    = var(parameter_ests(:));
end

% Inferred parameter values
Re = median(rmeans(1,ok));
Ra = median(rmeans(2,ok));
Rr = median(rmeans(3,ok));


fprintf('Model    excitation strength  is  %0.2e\n',model.rAA);
fprintf('Inferred excitation strength  is  %0.2e\n',Re);
fprintf('Model    A->R transition rate is  %0.2e\n',model.linearRates(2));
fprintf('Inferred A->R transition rate is  %0.2e\n',Ra);
fprintf('Model    R->Q transition rate is  %0.2e\n',model.linearRates(3));
fprintf('Inferred R->Q transition rate is  %0.2e\n',Rr);

learningRate = 0.8;
model.rAA                = learningRate.*Re+(1-learningRate)*model.rAA;
model.linearRates(2)     = learningRate.*Ra+(1-learningRate)*model.linearRates(2);
model.linearRates(3:end) = learningRate/2*Rr...
                         + (1-learningRate/2)*model.linearRates(3:5);

exc   = cell2mat(infe);
maxRe = 0.5./(dt*prctile(exc(:),90))
mode.rAA = min(maxRe,model.rAA);



