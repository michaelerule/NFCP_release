function [solution, fit, emuFit, info] = GPUCB_gpml(objectiveFunc, dim, searchLBs, searchUBs, options)
    %{
    GPUCB Perform optimisation using the GP UCB algorithm
    [solution, fit, emuFit, info] = GPUCB_gpml(objectiveFunc, dim, searchLBs, searchUBs, options)

    GPUCB implementation that makes use of the RW book implementation
    can be downloaded at:
    http://www.gaussianprocess.org/gpml/code/matlab/doc/

    Dimitrios Milios (dmilios@inf.ed.ac.uk)

    Parameters
    ----------
    objectiveFunc: handler to the objective function
    dim:           dimension of the search-space
    searchLBs:     vector containing the lower bounds of the search
    searchUBs:     vector containing the upper bounds of the search
    options:       struct with the optimisation options

    Returns
    -------
    solution: vector containing the best solution found
    fit:      the actual fitness of the given solution
    emuFit:   the emulated fitness (by the GP) of the given solution
    into:     struct that contains information about the search
    %}

	info.amplitude = 0;
	info.lengthscale = 0;
	info.totalEvaluations = 0;
	info.elapsedHyperparams = 0;
	info.elapsedGPUCB = 0;

	if options.sparseGP,
		if options.initSampleNumber<options.inducingPoints,
		    warning('initSampleNumber cannot be less than inducingPoints for sparseGP');
		    options.inducingPoints = options.initSampleNumber;
	    end
    end

	initSampleNumber = options.initSampleNumber;
	gridSampleNumber = options.gridSampleNumber;
	maxIterations    = options.maxIterations;
	maxEvaluations   = options.maxEvaluations;
	sigma2           = options.sigma2;
	beta             = options.beta;
	
    % Support user-supplied pre-initialization of data points
    if isfield(options,'presampled') && options.presampled,
        if ~isfield(options,'userSuppliedSamplesX'),
            error('If using the presampled option, please also supply the presampled point locations in options.userSuppliedSamplesX');
        end
        if ~isfield(options,'userSuppliedSamplesY'),
            error('If using the presampled option, please also supply the presampled point locations in options.userSuppliedSamplesY');
        end
        userSuppliedSamplesX = options.userSuppliedSamplesX;
        userSuppliedSamplesY = options.userSuppliedSamplesY(:);
        preSampleNumber = numel(userSuppliedSamplesY);
        if ~preSampleNumber==size(userSuppliedSamplesX,1),
            error('userSuppliedSamplesX does not match userSuppliedSamplesY in size, or it is transposed');
        end
    else
        preSampleNumber = 0;
    end
        
	% preallocate memory for training set -------------------------
	preallocmemX = zeros(maxEvaluations + initSampleNumber + preSampleNumber, dim);
	preallocmemY = zeros(maxEvaluations + initSampleNumber + preSampleNumber, 1);
	soFar = initSampleNumber + preSampleNumber;
	% -------------------------------------------------------------

    % Assign user-supplied pre-initialization
    if preSampleNumber>0,
        preallocmemX(1:preSampleNumber, :) = userSuppliedSamplesX;
        preallocmemY(1:preSampleNumber)    = userSuppliedSamplesY;
    end

	% initial grid sampling ---------------------------------------
	searchRanges = searchUBs - searchLBs;
	rangesMatrix = repmat(searchRanges, initSampleNumber, 1);
	lowerBoundsMatrix = repmat(searchLBs, initSampleNumber, 1);
	preallocmemX(preSampleNumber+1:soFar, :) = ...
	    rand(initSampleNumber, dim) .* rangesMatrix + lowerBoundsMatrix;
	for i = 1:initSampleNumber
		preallocmemY(preSampleNumber+i) = objectiveFunc(preallocmemX(i, :));
	end
	% Copy portion of points used thusfar
	gpX = preallocmemX(1:soFar, :);
	gpY = preallocmemY(1:soFar);

	rangesMatrix = repmat(searchRanges, gridSampleNumber, 1);
	lowerBoundsMatrix = repmat(searchLBs, gridSampleNumber, 1);

	
	% initialize GP ----------------------------------------------
	meanfunc = @meanConst; hyp.mean = 0;
	ell = ones(1, dim); sf = 1.0; hyp.cov = log([ell sf]);
	hyp.lik  = log(sigma2);
	likfunc = @likGauss;
	if ~options.sparseGP
		covfunc = @covSEard;
		inferenceFunc = @infExact;
	else
		nu = options.inducingPoints;
		inducingPoints = rand(nu, dim) .* repmat(searchRanges, nu, 1) + repmat(searchLBs, nu, 1);
		covfunc = {@covFITC, {@covSEard}, inducingPoints};
		inferenceFunc = @infFITC;
	end
	
	% hyperparam optimisation
	tic; % ====================================================================
	hyp = minimize(hyp, @gp, -40, inferenceFunc, meanfunc, covfunc, likfunc, gpX, gpY);
	info.elapsedHyperparams = toc;  % =========================================
	disp(['Hyperparam opt: ' num2str(info.elapsedHyperparams) ' sec'])
	
	info.amplitude = exp(hyp.cov(end));
	info.lengthscale = exp(hyp.cov(1:dim));
	info.sigma2 = log(sqrt(sigma2));

	% GPUCB 
	tic; % ====================================================================
	
	[~, ~, mu, var, ~] = gp(hyp, inferenceFunc, meanfunc, covfunc, likfunc, gpX, gpY, gpX);
	decisionObs = mu + beta .* sqrt(var);

	localOptimOptions = optimset('Display', 'off', 'Algorithm', 'active-set');
	localOptArea = 0.1;

	iteration = 0;
	evaluation = 0;
	msg_n = 0;
	while iteration < maxIterations && evaluation < maxEvaluations
	
		gpTest = rand(gridSampleNumber, dim) .* rangesMatrix + lowerBoundsMatrix;
		[~, ~, mu, var, ~] = gp(hyp, inferenceFunc, meanfunc, covfunc, likfunc, gpX, gpY, gpTest);
		
		decision = mu + beta .* sqrt(var);
		
		[maxDec, maxDecIdx] = max(decision);
		candidate = gpTest(maxDecIdx, :);
		
		localLBs = candidate - localOptArea * searchRanges / 2;
		localUBs = candidate + localOptArea * searchRanges / 2;
		gpUBFn = @(x) gpUpperBound_neg(x, gpX, gpY, hyp, beta);
		[candidate, cost] = fmincon(gpUBFn, candidate, [], [], [], [],...
			localLBs, localUBs, [], localOptimOptions);
		maxDec = -cost;
		
		if maxDec > max(decisionObs)
			newInstance = candidate;
			newObs = objectiveFunc(newInstance);
			soFar = soFar + 1;
			preallocmemX(soFar, :) = newInstance;
			preallocmemY(soFar) = newObs;
			gpX = preallocmemX(1:soFar, :);
			gpY = preallocmemY(1:soFar);

			[~, ~, mu, var, ~] = gp(hyp, inferenceFunc, meanfunc, covfunc, likfunc, gpX, gpY, gpX);
			
			decisionObs = mu + beta .* sqrt(var);
			evaluation = evaluation + 1;
		end
		
		iteration = iteration + 1;
		msg = ['iteration: ' num2str(iteration) ',  evaluations: ' num2str(evaluation)];
		fprintf(repmat('\b', 1, msg_n));
		fprintf(msg);
		msg_n = numel(msg);
	end

	[~, ~, mu, ~, ~] = gp(hyp, inferenceFunc, meanfunc, covfunc, likfunc, gpX, gpY, gpX);
	info.elapsedGPUCB = toc;  % =========================================================
	
	[emuFit, solutionIndex] = max(mu);
	solution = gpX(solutionIndex, :);
	fit = gpY(solutionIndex);
	
	localLBs = solution - localOptArea * searchRanges / 2;
	localUBs = solution + localOptArea * searchRanges / 2;
	gpMuN = @(x) gpMean_neg(x, gpX, gpY, hyp);
	[solution, cost] = fmincon(gpMuN, solution, [], [], [], [],...
		localLBs, localUBs, [], localOptimOptions);
	fit = -cost;
	
	info.totalEvaluations = evaluation + initSampleNumber;
	
	
	function mu_neg = gpMean_neg(x, gpX, gpY, hyp)
		mu_neg = -gp(hyp, inferenceFunc, meanfunc, covfunc, likfunc, gpX, gpY, x);
	end
	
	function up_neg = gpUpperBound_neg(x, gpX, gpY, hyp, beta)
		[mu, var] = gp(hyp, inferenceFunc, meanfunc, covfunc, likfunc, gpX, gpY, x);
		up_neg = -(mu + beta .* sqrt(var));
	end
	
end
