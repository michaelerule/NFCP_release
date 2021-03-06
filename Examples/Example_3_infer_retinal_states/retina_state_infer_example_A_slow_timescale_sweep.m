%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up workspace
close all; clear all;
addpath('./','../','../../','../../NFCP'); 
NFCP_init

% Set this to the path of NFCP_datasets on your system
% (Contact mrule7404@gmail.com for a copy of NFCP_datasets)
% Here is a dropbox share link (current as of 2019/08/18)
% https://www.dropbox.com/sh/apifkoowouvq2jc/AABLysnAJWs1YY2fmB4HKVRWa?dl=0
datadir  = '/home/mer49/Dropbox (Cambridge University)/Datasets/NFCP_datasets/';

% Script options
SRES      = 10; % spatial resolution (10, 15 and 20 are available)
df        = 'P11_binned_100ms_20150219_leptin_retina1_control';

% Locate data file
subdir   = 'binned_spikes/spatiotemporal_binned/renewal_model/version_4_selected_datasets/';
fn       = sprintf('%d/%s_resolution_%d',SRES,df,SRES);
suffix   = '_preinitialized_using_renewal_model.mat';
filepath = [datadir subdir fn suffix];
loadarchive 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model configuration

% Heuristic model configuration tuned for P11_binned_100ms_20150219_leptin_retina1_control

model.rAA           = 12;    % Excito-Excitatory interaction strength
model.thr           = 0.12;  % Activation threshold
model.sigma         = 0.1;   % Standard-deviation for excitatory interaction kernel
model.description = [ ...
%    Q  A  R   rate
    -1  1  0   0.0 % spontaneous excitation 
     0 -1  1   1.8 % enter slow refractory loop
     1  0 -1   0.1 % slow refractory recovery
     ];

% Algorithim configuration parameters
model.dt            = 0.1;   % time interval between observations
model.n             = SRES;  % Simulation grid size
model.cutoff        = true;  % Remove scales finer than interaction radius?
model.ini_state_var = 1e-0;  % Initial variance in state estimate
model.reg_state_var = 1e-2;  % Added uncertainty in state 
model.ini_count_var = 0;     % Initial variance for neuron count
model.reg_count_var = 0;     % Added uncertainty in neuron count
model.reg_inverse   = 1e-12; % Diagonal regularizer for matrix inversions
model.ini_reg_diag  = 1e-12; % Initial variance added to diagonal
model.reg_diag      = 1e-12; % Added diagonal regularizaer
model.reg_precision = 1e-12; % Regularization for precision matrix
model.update        = 'Laplace'; % Measurement update method
model.minrate       = 1e-6;  % Regularization parameter to avoid zeros
model.maxiter       = 1000;  % Limit Newton-raphson measurement iterations
model.tol           = 1e-6;  % Convergence tolerance for Newton-Raphson
model.verbosity     = 0;     % How many logging details to print
model.safety        = 0;     % Whether to do extra numeric stability checks
model.alpha         = 1;     % Dispersion paramter, 1=Poisson
effectivepopsize    = 100;   % Effective population size for scaling noise

% Colors and rescaling for plotting
model.cscale        = [1 5 1];
model.colors        = [0 0 1;  % Q: blue
                       1 0 0;  % A: red
                       0 1 0]; % R: green
model.names = strsplit('Q A R');
ratescale = model.cscale(2);
model.names{2} = sprintf('A×%d',ratescale);

model.ss_rescale    = 1/effectivepopsize;
fprintf('model.ss_rescale = %f\n',model.ss_rescale);
model = initializeModel(model);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Re-compute gains/biases
% Create a histogram of 2D xypoints (the observation counts)
% These are premultiplied by alpha
counts = {};
for i=1:size(xydata,2), counts{i} = binCounts(model,xydata{i}); end
counts = cell2mat(counts);
for i=1:N*N, b2(i) = prctile(counts(i,1:end),50); end
for i=1:N*N, g2(i) = prctile(counts(i,1:end),99.999); end
goodchannels = inhomogeneous_gain>0;
g2 = g2(:);
model.gain = mean(g2(goodchannels)./inhomogeneous_gain(goodchannels))/model.volume;
model.bias = b2;
model.inhomogeneous_gain = gaussian2DblurOperator(model.n,model.sigma.*model.n)*model.inhomogeneous_gain;
model.bias = gaussian2DblurOperator(model.n,model.sigma.*model.n)*model.bias';
fprintf('Old gain was %0.4f\n',gain);
fprintf('New gain is  %0.4f\n',model.gain);
fprintf('Old bias was %0.4f\n',mean(bias(goodchannels)));
fprintf('New bias is  %0.4f\n',mean(model.bias(goodchannels)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Soft boundary around active regions
% (messy code gets reflected boundary conditions for a spatial blur)
for i=1:N*N, m2(i) = mean(counts(i,1:end)); end
goodregions = m2>prctile(m2,15);
U = reshape(goodregions,N,N);
Upadded = [U fliplr(U)];
Upadded = [Upadded; flipud(Upadded)];
K       = floor(N/2);
Upadded = circshift(Upadded,[K,K]);
Upadded = gaussian2DblurOperator(N*2,0.1*N)*Upadded(:);
Upadded = reshape(Upadded,N*2,N*2);
U       = Upadded(K:K+N-1,K:K+N-1);
U       = U(:);
ini = [U*0.7 U*0 U*0.3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordinate sweep the slow R->Q rate

START = 1900;
STOP  = 4000;
model.update     = 'Laplace'; % Measurement update method
model.likemethod = 'ELBO'; K=1;

figure(6); clf;


Nscan  = 4*12+1;
rrates = 10.^linspace(-3,1,Nscan);
obj    = @(m) modelLikelihood(m,ini,xydata(:,START:STOP));
par2m  = @(r) initializeModel(applyOptions(model,{'linearRates',[0 model.linearRates(2) r]},false));
all_lr = {};
parfor i=1:Nscan, all_lr{i} = obj(par2m(rrates(i))); end
all_lr = reshape(cell2mat(all_lr),K,Nscan);
save('RQ_ll.mat','all_lr');

subplot(221); cla;
llb2 = log2(exp(all_lr));%.*rrates.*log(10));
semilogx(rrates,llb2-max(llb2)); hold on;
fprintf(1,'Optimal rate is %f\n',rrates(argmax(all_lr)));
[llmax,iimax] = max(all_lr);
best = rrates(iimax); 
xline(rrates(iimax));
fprintf(1,'Optimal rate is %f\n',best);
xlabel('R→Q rate (transitions/s)');
ylabel('log_2(likelihood) + const. (bits/sample)');
ylim([-7,1])



Nscan  = 4*10+1;
qrates = 10.^linspace(-2,2,Nscan);
obj    = @(m) modelLikelihood(m,ini,xydata(:,START:STOP));
par2m  = @(q) initializeModel(applyOptions(model,{'linearRates',[0 q model.linearRates(3)]},false));
all_lq = {};
parfor i=1:Nscan, all_lq{i} = obj(par2m(qrates(i))); end
all_lq = reshape(cell2mat(all_lq),K,Nscan);
save('AR_ll.mat','all_lq');

subplot(222); cla;
llb2 = log2(exp(all_lq));%.*rrates.*log(10));
semilogx(qrates,llb2-max(llb2)); hold on;
best = qrates(argmax(all_lq));
xline(best);
fprintf(1,'Optimal rate is %f\n',best);
xlabel('A→R rate (transitions/s)');
ylabel('log_2(likelihood) + const. (bits/sample)');
ylim([-7,1])



Nscan  = 4*2+1;
arates = 10.^linspace(-2,2,Nscan);
obj    = @(m) modelLikelihood(m,ini,xydata(:,START:STOP));
par2m  = @(a) initializeModel(applyOptions(model,{'rAA',a},false));
all_la = {};
parfor i=1:Nscan, all_la{i} = obj(par2m(arates(i))); end
all_la = reshape(cell2mat(all_la),K,Nscan);
save('QA_ll.mat','all_la');

subplot(223); cla;
llb2 = log2(exp(all_la));%.*rrates.*log(10));
semilogx(arates,llb2-max(llb2)); hold on;
best = arates(argmax(all_la));
xline(best);
fprintf(1,'Optimal rate is %f\n',best);
xlabel('QA→AA rate (transitions/s/cell)');
ylabel('log_2(likelihood) + const. (bits/sample)');
ylim([-7,1])





