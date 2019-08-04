% Set up workspace
close all; clear all; rng('shuffle');      % Clean up namespace
addpath('./','../','../../NFCP'); % Add library path (relative path)
NFCP_init

% Need to test EP update with various link functions, some noise, etc. 
n    = 17;
nstates = 3;
m    = n*n;
dt   = 1.0;
dx   = 1/m;
gain = (30+rand(m,1)*0.001).*dx.*dt;
%gain = ones(m,1).*dx.*dt*30;
bias = rand(m,1)*0.001;
%bias = zeros(m,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct minimal model struct needed to get us off the ground
model.n  = n;
model.nn = m;
model.dt = dt;
model.dx = dx;
model.safety = 2;
model.premultiplied_gain = gain;
model.adjusted_bias = bias;
model.link = 'squared';
model.llreg = 1e-9;
% We don't really use the dynamics, but define them so initializeModel is happy
model.rQA     = 2e-1;  % Rate of spontaneous excitation (modeled as shot noise)
density       = 50;    % Number of "neurons" per unit area
model.rAA     = 1.40./density;  % Excito-Excitatory interaction strength
model.thr     = 8e-3;  % Nonzero threshold for depolarization
model.sigma   = 0.1;   % Standard-deviation for excitatory interaction kernel
model.alpha   = 1;     % Dispersion paramter, 1=poisson
model.description = [ ...
%    Q  A  R1   rate
    -1  1   0   model.rQA     % spontaneous excitation (zero because handled by rQA)
     0 -1   1   4e-1  % slow refractory loop
     1  0  -1   32e-4 % slow refractory recovery
     ];
model = initializeModel(model);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample something with a little spatiotemporal structure

% Random intensity fields with some spatial correlations
M    = rand(m*3,1);    % Random prior
A    = randn(m*3)*0.1; % Prior covariance
A    = model.K2D3*A;    % Add spatial correlation to covariance factor
M    = model.K2D3*M;    % Add spatial correlation to mean
C    = A*A';     

% Random state
x    = max(mvnrnd(M,C)',0);

% Generate an incorrect prior for test
Mq   = rand(m*3,1);     % Random prior
Aq   = randn(m*3);      % Prior covariance
Aq   = model.K2D3*A;    % Add spatial correlation to covariance factor
Mq   = model.K2D3*M;    % Add spatial correlation to mean
Cq   = A*A';    

% Active subspace
xa   = x(m+1:2*m);
Ma   = Mq(m+1:2*m);
Ca   = Cq(m+1:2*m,m+1:2*m);
Pa   = cinv(Ca);

% Poisson sample, linear link
rate = gain.*xa.^2+bias;
y = poissrnd(rate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check that likelihood surfaces look OK

b  = max(Ma+10*diag(Ca).^0.5);
xp = linspace(1e-9,b,400);
subplot(332)
model.link = 'log';
logpyx = likelihoodSurface(model,xp,y);
plot(logpyx')
subplot(331)
model.link = 'linear';
logpyx = likelihoodSurface(model,xp,y);
plot(logpyx')
subplot(333)
model.link = 'squared';
logpyx = likelihoodSurface(model,xp,y);
plot(logpyx')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check that prior surfaces look OK

V = diag(C);

subplot(334)
logpx = priorSurface(M,1./V,xp);
plot(exp(logpx'));

% Check that truncated Gaussian is matching moments accurately
p = normrows(exp(logpx));        
% Integrate to get posterior moments
trM = sum(bsxfun(@times,p',xp')',2);
trV = sum(bsxfun(@times,p,(outer(xp,trM,@minus)).^2'),2);

rms(trM(:)-M(:))
rms(trV(:)-V(:))

% This appears to work? 
Vadj = V/(1-2/pi);
Madj = M - sqrt(Vadj)*sqrt(2/pi);
samples = bsxfun(@plus,bsxfun(@times,randn(867,10000),sqrt(Vadj)),Madj);
samples(samples<0) = NaN;
samples(samples>b) = NaN;
v2 = nanvar(samples,[],2);
m2 = nanmean(samples,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check that posterior surfaces look OK

[pxy,pM,pV,pP,dM,dP,dPM] = integrateMarginals(Ma,diag(Ca),xp,logpyx);
subplot(335)
plot(pxy');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[loglikelihood,M2,C2,Ma2,Pa2,rMa,rPa,rCa] = measurementEP(model,Mq,Cq,y);

subplot(336)
rate0 = Ma(:).^2.*gain + bias;
scatter(rate,rate0)
subplot(337)
rate2 = Ma2(:).^2.*gain + bias;
scatter(rate,rate2)
subplot(338)
logpx = priorSurface(Ma2,1./diag(cinv(Pa2)),xp);
plot(exp(logpx'));
subplot(339)
scatter(xa,Ma2(:))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image plots

figure(2);
subplot(241);
imagesc(reshape(xa,n,n))
title('True');
subplot(242);
imagesc(reshape(Ma,n,n))
title('Incorrect Prior');
subplot(243);
imagesc(reshape(Ma2,n,n))
title('EP');
% Compare to Laplace
[ll,M3,C3] = measurementSubspace(model,M,C,y);
Ma3 = M3(m+1:2*m);
subplot(244);
imagesc(reshape(Ma3,n,n));
title('Laplace');

% Compare, cell by cell, the reconsructed means
subplot(245)
scatter(Ma2,Ma3);
b = max(max(Ma2),max(Ma3));
xlim([0 b]);
ylim([0 b]);
xlabel('EP');
ylabel('Laplace');
subplot(246)
scatter(xa,Ma);
xlabel('True');
ylabel('Prior');
subplot(247)
scatter(xa,Ma2);
xlabel('True');
ylabel('EP');
subplot(248)
scatter(xa,Ma3);
xlabel('True');
ylabel('Laplace');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare to effects of truncation

figure(3)

% Truncation-only update
[loglikelihood,M4,C4,Ma4,Pa4,rPa4,rCa4] = measurementEP(model,Mq,Cq,y.*NaN);

subplot(241);
imagesc(reshape(Ma4,n,n))
subplot(245);
scatter(Ma,Ma4);
xlabel('Prior');
ylabel('Truncation');
upper = max(max(Ma),max(Ma4));
xlim([0 upper]);
ylim([0 upper]);


% Manually extract the truncation effect
dPt  = Pa4 - Pa;
dPMt = Pa4*Ma4 - Pa*Ma;
dMt  = dPt\dPMt;

% Manually remove the truncation effect from the EP posteroir
cPa  = Pa2 - dPt;
cPMa = Pa2*Ma2 - dPt*dMt;
cMa  = cPa\cPMa;

subplot(242);
imagesc(reshape(cMa,n,n))
subplot(246);
scatter(cMa,Ma3);
xlabel('Corrected');
ylabel('Laplace');
upper = max(max(cMa),max(Ma3));
xlim([0 upper]);
ylim([0 upper]);

subplot(243);
imagesc(reshape(cMa,n,n))
subplot(247);
scatter(xa,cMa);
xlabel('True');
ylabel('Corrected');
upper = max(max(cMa),max(xa));
xlim([0 upper]);
ylim([0 upper]);

subplot(244);
imagesc(reshape(dMt,n,n))
subplot(248);
scatter(xa,dMt);
xlabel('True');
ylabel('Correction');




figure(5)

% This appears to work? 
Vadj = V/(1-2/pi);
Madj = M - sqrt(Vadj)*sqrt(2/pi);
samples = bsxfun(@plus,bsxfun(@times,randn(867,10000),sqrt(Vadj)),Madj);
samples(samples<0) = NaN;
samples(samples>b) = NaN;
v2 = nanvar(samples,[],2);
m2 = nanmean(samples,2);

subplot(221)
scatter(M,trM);
subplot(222)
scatter(V,trV);
subplot(223)
scatter(M,m2);
subplot(224)
scatter(V,v2);


