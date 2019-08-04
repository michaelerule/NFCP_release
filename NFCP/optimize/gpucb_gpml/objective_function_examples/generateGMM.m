function [weights, means, Sigmas] = generateGMM(k, d, lbound, ubound, seed)
% generateGMM Randomly generate a Gaussian Mixture Model
% [weights, means, Sigmas] = generateGMM(k, d, lbound, ubound, seed)
%
%	k:      the number of components for the GMM
% 	d:      the dimension of the space (i.e. d-variate GMM)
%	lbound: vector containing the lower bounds
%	ubound: vector containing the upper bounds
%	seed:   the seed to be used by the random engine
%
%	weights: the component weights for the GMM
%	means:   the component means for the GMM
%	Sigmas:  the component (diagonal) covariance matrices for the GMM
%
% Dimitrios Milios (dmilios@inf.ed.ac.uk)

	if ~isvector(lbound) || ~isvector(ubound)
		error('lbound and ubound shoud be vectors');
	end

	if length(lbound) == 1 && length(ubound) == 1
		range = ubound - lbound;
		lbounds = ones(1, d) .* lbound;
		ubounds = ones(1, d) .* ubound;
		ranges = ubounds - lbounds;
	else 
		if length(lbound) == d && length(ubound) == d
			lbounds = lbound;
			ubounds = ubound;
			ranges = ubounds - lbounds;
		else
			error('lbound and/or ubound are incompatible with the dimension d');
		end
	end

	% volume of d-ball:
	% v = range^d / 10;
	% r = v^(1/d) * gamma(d/2 + 1)^(1/d) / sqrt(pi);
	
	logv = sum(log(ranges)) - log(10);
	logr = (1/d) * logv + (1/d) * gammaln(d/2 + 1) - log(sqrt(pi));
	r = exp(logr);
	
	previousSeed1 = rand('seed');  % save previous state of the random generator
	rand('seed', seed);
	weights = rand(1, k);
	weights = weights / sum(weights);
	means = rand(k, d) .* repmat(ranges, k, 1) + repmat(lbounds, k, 1);
	Sigmas = zeros(d, d, k);
	for i = 1:k
		Sigmas(:, :, i) = diag(rand(1, d) .* r);
	end
	rand('seed', previousSeed1);  % load previous state of the random generator
end
