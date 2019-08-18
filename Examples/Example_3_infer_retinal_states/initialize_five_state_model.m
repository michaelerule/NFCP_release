% Excitation model options
model.rQA     = 0;    % Rate of spontaneous excitation (modeled as shot noise)
model.thr     = 1e-4;    % Threshold for depolarization
model.sigma   = 0.1;  % Standard-deviation for excitatory interaction kernel

% initialize_five_state_model.m
% Slow refractory model
model.names = strsplit('Q A R1 R2 R3');
model.description = [...
    -1  1  0  0  0 model.rQA % spontaneous excitation (zero because handled by rQA)
     0 -1  1  0  0 rA        % slow refractory loop
     0  0 -1  1  0 rR
     0  0  0 -1  1 rR
     1  0  0  0 -1 rR];    % slow refractory recovery

% Plotting colors as R G B triplets in range 0 to 1
model.colors = [0 0 1;1 0 0;0.25 1 0;0 1 0.25;0.125 1 0.125];

% Plotting configuration: color scales for mapping Q/A/R fields to display
A_color_scale = 4;
Q_color_scale = 1;
R_color_scale = 1;
model.cscale = [Q_color_scale A_color_scale R_color_scale R_color_scale];

% initial conditions of intensity fields:
Q = ones(model.n,model.n)*0.6; % Quiescent 
A = ones(model.n,model.n)*0.1; % Active
R1= ones(model.n,model.n)*0.1; % Recovery state
R2= ones(model.n,model.n)*0.1; % Recovery state
R3= ones(model.n,model.n)*0.1; % Recovery state
ini= [Q A R1 R2 R3];

% Total population normalized to 1, all rates refer to fractions
model.normalized    = true;
model.verbosity     = 0;  % do not print debugging information
model.safety        = 0;  % no not perform safety checks

% Measurement update and likelihood methods
model.update        = 'Laplace'; % Measurement update method
model.update        = 'Laplace-subspace'; % Measurement update method
model.likemethod    = 'ML'; % maximum likelihood approach discards the prior

% Hyperparameters: prevent singular matrices
model.minrate       = 1e-6;
model.cutoff        = false; % Remove scales finer than interaction radius?
model.reg_state_var = 1e-2;  % Added uncertainty in state 
model.reg_count_var = 1e-8;  % Added uncertainty in neuron count
model.reg_diag      = 1e-8;  % Added diagonal regularizaer
model.reg_inverse   = 1e-8;  % Diagonal regularizer for matrix inversions
model.ini_state_var = 1e-4 ; % Initial variance in state estimate
model.ini_count_var = 1e-8;  % Initial variance for neuron count
model.ini_reg_diag  = 1e-8;  % Initial variance added to diagonal
model.reg_precision = 1e-8;  % Regularization for precision matrix
model.normalized    = false;

model = initializeModel(model);
