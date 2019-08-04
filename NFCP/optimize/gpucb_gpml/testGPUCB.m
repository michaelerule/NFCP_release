% Example of usage for GPUCB
% Examples of objective functions can be found in 'objective_function_examples'
%
% Dimitrios Milios (dmilios@inf.ed.ac.uk)

addpath('gpml-matlab-v3.4-2013-11-11'); startup;
addpath('objective_function_examples');

%% Set the dimension of the search space and the corresponding bounds
d = 2;
searchLBs = ones(1, d) .* 0.001;
searchUBs = ones(1, d) .* 10;
searchRanges = searchUBs - searchLBs;

%% Create and objective function (i.e. GMM in this example)
% -----------------------------------------------
k = 5;
modelSeed = 12367;
[weights, means, Sigmas] = generateGMM(k, d, searchLBs, searchUBs, modelSeed);
objectiveFunc = @(x) gmmLogPDF(x, weights, means, Sigmas);  % gmmPDF  or  gmmLogPDF

% Find the true maximum
values = zeros(1, k);
for i = 1:k
	values(i) = objectiveFunc(means(i,:));
end
[~, maxIndex] = max(values);
trueMaximum = means(maxIndex, :);
trueMaximumFit = objectiveFunc(trueMaximum);
clear maxIndex;
% -----------------------------------------------

%% Set search options and perform the search via GPUCB
% -----------------------------------------------
methodSeed = cputime; % cputime
rand('seed', methodSeed);
randn('seed', methodSeed + 1);


options.initSampleNumber = 120;
options.latinHypercubeSampling = true;
options.latinHypercube_subdivisions = 5;
options.latinHypercube_samplesPerSubspace = 4;

options.sparseGP = true;
options.inducingPoints = 20; % number of inducing points for sparse GP (it should be smaller that the training set size)

options.gridSampleNumber = 100;
options.maxIterations = 100;
options.maxEvaluations = 100;
options.maxAddedPointsNoImprovement = 100;
options.improvementFactor = 1.01;
options.maxFailedAttempts = 100;
options.debugPlots = 1;

options.sigma2 = 0.05;  % noise parameter for the GP
options.beta = 2;       % beta parameter for the GP-UCB (i.e. upper quantile = beta * stdev)

[solution, fit, emuFit, info] = GPUCB_gpml(objectiveFunc, d, searchLBs, searchUBs, options);



%% Show results
disp(' ');
disp(['  amplitude used: ' num2str(info.amplitude)])
disp(['lengthscale used: ' num2str(info.lengthscale)])
disp(['     sigma2 used: ' num2str(info.sigma2)])
disp('------------------------------------------------')
disp(['Best solution fitness: ' num2str(fit)])
disp(['     Emulated fitness: ' num2str(emuFit)])
disp('------------------------------------------------')
disp(['Distance from actual best: ' num2str(norm(solution - trueMaximum))])
disp('------------------------------------------------')
disp(['  Actual best fitness: ' num2str(trueMaximumFit)])
disp('------------------------------------------------')
disp(['Evaluations: ' num2str(info.totalEvaluations)])
disp(['Hyperparam opt: ' num2str(info.elapsedHyperparams) ' sec'])
disp(['GPUCB:          ' num2str(info.elapsedGPUCB) ' sec'])


