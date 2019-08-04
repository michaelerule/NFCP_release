#!/usr/bin/env python
# -*- coding: UTF-8 -*-
def initializeModel(model,varargin):
    r'''

    Precompute commonly used values and initialize the model.

    +--------------------------+-------------------------------------------+
    | Variable                 | Units                                     |
    +==========================+===========================================+
    | model.dt                 | seconds      /   bin_t                    |
    +--------------------------+-------------------------------------------+
    | model.n                  | bin_x        /   array_x                  |
    +--------------------------+-------------------------------------------+
    | model.nn                 | bin_x²       /   array_x²                 |
    +--------------------------+-------------------------------------------+
    | model.dx                 | array_x²     /   bin_x²                   |
    +--------------------------+-------------------------------------------+
    | model.volume             | s∙array_x²   /   bin_x²∙bin_t             |
    +--------------------------+-------------------------------------------+
    | model.bias               | spikes       /   bin_x²∙bin_t             |
    +--------------------------+-------------------------------------------+
    | model.gain               | spikes       /   s∙array_x²               |
    +--------------------------+-------------------------------------------+
    | model.alpha              | spike mean   /   spike standard dev.      |
    +--------------------------+-------------------------------------------+
    | model.sigma              | array_x (1=linear size of array)          |
    +--------------------------+-------------------------------------------+
    | model.linearRates        | 1/neuron/s                                |
    +--------------------------+-------------------------------------------+
    | model.rAA                | 1/neuron²/s (1/fraction²/s if normalized) |
    +--------------------------+-------------------------------------------+
    | model.inhomogeneous_gain | unitless                                  |
    +--------------------------+-------------------------------------------+
    | model.adjusted_bias      | unitless                                  |
    +--------------------------+-------------------------------------------+
    | model.adjusted_gain      | alpha∙spikes /   bin_x²∙bin_t             |
    +--------------------------+-------------------------------------------+
    | rate (inside obj. fun.)  | spikes       /   bin_x²∙bin_t             |
    +--------------------------+-------------------------------------------+

    Example
    -------
    Please see example scripts for details; this is an example of the
    a defined linear transition model. ::

        % Transition rates
        model.names = strsplit('Q A R1 R2');

        %    Q  A  R1  R2  rate
        model.description = [
            -1  1   0   0  2e-1 % spontaneous excitation
             0 -1   1   0  3e-1 % slow refractory loop
             0 -1   0   1  1e-1 % fast refractory loop
             1  0  -1   0  3e-4 % slow refractory recovery
             1  0   0  -1  1e-2 % fast refractory recovery
             ];

    Parameters
    ----------
    model : struct

    Other Parameters
    ----------------
    method : `string`, default 'momentClosure'
        State inference method to use. Can be `'momentClosure'` or `'LNA'`
    dt : `float`, default 1.0
        Time resolution of sampled; For state or model inference, point-
        process data should be binned at resolution `dt`
    n : `int`, default 10
        Resolution of the `n x n` grid of basis functions onto which the
        continuous system is projects. Recommended no smaller than `7 x 7`
    gain : `float`, default 1
        Global model gain parameter. By convention the second model state
        is interpreted as the active (A) state. The intensity of this state
        is multiplied by `gain` to estimate the expected number of events
        (spikes) observed.
    bias : `float` or `vector`, default 0
        Per-channel region in the point-process intensity. Can be a scalar.
        If a vector, should have one entry for every region in the `n x n`
        simulation grid.
    alpha : `float`, default 1
        Parameter controlling the dispersion of the count process.
        `alpha=1` is Poisson. Equal to the reciprocal of the fano factor.
    cutoff : `bool`, default true
        Whether to let the interaction length (`model.sigma`) set the
        minimum spatial resolution of the model. If `true`, spatial
        frequencies will be attenuated by `model.sigma`.
    reg_state_var : `float`, default 1e-3
        Covariance matrix regularization parameter (additive noise).
        Regularization orthogonal to neuron counts and spatial correlations
    reg_count_var : `float`, default 1e-8
        Covariance matrix regularization parameter (additive noise).
        Regularization orthogonal to spatial correlations
    reg_diag : `float`, default 1e-8
        Covariance matrix regularization parameter (additive noise).
        Diagonal regularization'
    reg_inverse : `float`, default 1e-3
        Precision matrix regularization. Should remain small to avoid
        under-estimating model variance.
        Diagonal regularization for recovering covariances for moment closure
    interpolate : `bool`, default true
        Whether to use interpolation when binning spikes to basis elements
    ss_rescale : `float`, default 1.0
        Effective system-size re-scaling; multiplies noise levels. Combine
        with `model.normalized=true` to run inference using normalized
        population densities, with the noise adjusted by `ss_rescale`
    maxiter : `int`, default 20
        Newton-Raphson iteractions in measurement update
    link : `string`, default 'linear'
        Reserved option for incorporating new observation models;
        onlt `linear` is supported at the moment.
    llreg : `float`, default 0.0
        Default is no regularization for log-Determinant in likelihood;
        These should be positive-definite anyway.
    tol : `float`, default 1e-6
        Measurement convergence tolerance
    minrate : `float`, default 1e-6
        Minimum intensity in measurement update
    L2penalty : `float`, default 1e-14
        L2 penalty on intensities of model
    sqrtform : `bool`, default true
        Whether to use square-root filter updates when possible
    normalized : `bool`, default false
        Whether local intensities should be constrained to always sum to 1
    quantize : `bool`, default true
        Round parameters to finite sigfigs to facilitate caching
    nsigfig : `int>0`, default 3
        No sigfigs; smaller=more cache retrievals, but less accurate
    safety : `int`, default 1
        Safety level; 0=unchecked, 1=patch, 2=patch and warn, 3=abort
    verbosity : `int`, default 0
        Logging level; 0=nothing, 1=limited, 2=verbose, 3=very verbose
    ini_state_var : `float`, default `model.reg_state_var`
        Initial covariance, uncertainty in state
    ini_count_var : `float`, default `model.reg_count_var`
        Initial covariance, uncertainty in population count
    ini_reg_diag : `float`, default `model.reg_diag`
        Initial covariance, diagonal regularizer
    reg_precision : `float`, default `model.reg_diag`
        Setting precision regularization to match covariance reguarlization
        TODO: this doesn't sound right; FIX

    Returns
    -------
    struct
        Model struct with auxiliary variables initialized.

    '''
    pass#SKIPME
    '''#STARTCODE

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set defaults for optinal variables if they are not yet defined

    if numel(varargin)>0,
        model = applyOptions(model,varargin);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mandatory Parameters
    % Failure to specify model rate parameters, threshold, or interaction length is an error
    if ~isfield(model,'rAA'  ),
        error('Please define model.rAA, the strength of excito-excitatory interactions');
    end
    if ~isfield(model,'thr'  ),
        error('Please define model.thr, the finite excitation threshold');
    end
    if ~isfield(model,'sigma'),
        error('Please define model.sigma, the lateral interaction length');
    end

    %|=== PARAMETER ===|= DEFAULT VALUE =|===== NOTE ======|
        % Optional Parameters: State space method and spatiotemporal grid have defined defaults
        'method'       ,'momentClosure','(can also be LNA)'
        'dt'           ,1.0      ,''
        'n'            ,10       ,''
        % Optional Parameters: Observation model defaults to tracking density of Active cells
        'gain'         ,1        ,''
        'bias'         ,0        ,''
        'gamma'        ,1        ,'Exponentiate biases; used to tune log-link model'
        'alpha'        ,1        ,'(Poisson).'
        % Optional Parameters: regularization
        'cutoff'       ,false    ,'(interaction length sets minimum scale).'
        'reg_state_var',1e-5     ,'Regularize orthogonal to pop. size & spatial correlation'
        'reg_count_var',1e-8     ,'Regularization orthogonal to spatial correlations'
        'reg_diag'     ,1e-8     ,'Diagonal regularization'
        'reg_inverse'  ,1e-5     ,'Diagonal regularization for covariance inverse'
        'interpolate'  ,true     ,'Use interpolation when binning spikes to basis elements?'
        'ss_rescale'   ,1.0      ,'Inverse system-size parameter; multiplies noise levels'
        % Measurement update configuration
        'link'         ,'linear' ,'link function for observation'
        'update'       ,'Laplace','Measurement update algorithm; only Laplace supported presently'
        'likemethod'   ,'ML'     ,'only maximum-likelihood (ML) is supported at this time'
        % Newton-Raphson convergence control and regularization parameters
        'llreg'        ,0.0      ,'No regularization for log-Determinant in likelihood'
        'tol'          ,1e-4     ,'Measurement convergence tolerance'
        'minrate'      ,1e-4     ,'Minimum intensity in measurement update'
        'L2penalty'    ,1e-12    ,'L2 penalty on intensities of model (unused; TODO: remove)'
        % Other runtime parameters
        'sqrtform'     ,false    ,'Use square-root filter updates when possible?'
        'normalized'   ,false    ,'Constrain intensities to always sum to 1?'
        'safety'       ,1        ,'Safety 0=none, 1=patch, 2=warn, 3=abort'
        'verbosity'    ,0        ,'Logging 0=none, 1=some, 2=verbose, 3=highest'
        'cscale'       ,false    ,'Separate scales for each species'
        };
    model = applyOptions(defaults, model);
        % Optional Parameters: regularization
        'ini_state_var',model.reg_state_var,'Using regularization noise covariance as initial covariance'
        'ini_count_var',model.reg_count_var,'Using regularization noise covariance as initial covariance'
        'ini_reg_diag' ,model.reg_diag     ,'Using regularization noise covariance as initial covariance'
        'reg_precision',model.reg_diag     ,'Setting precision regularization to match covariance reguarlization'
        };
    model = applyOptions(defaults, model);

    % Check and warn if there are untested/unstable settings
    % Check the likelihood estimator
    % We've had some trouble with the prior contribution introducing
    % a significant bias in parameter estimates. We believe this is due
    % to the model being weakly identifiable, owing to sloppiness and a
    % large unobserved state space. The effect of the prior is to bias
    % parameters in a way that might not respect the true system
    if ~strcmp(model.likemethod,'ML'),
        warning(['Only maximum likelihood (ML) supported; "' model.likemethod '" untested']);
    end

    % Logging levels :
    if     model.safety<=1, model.error = @(s)[];   % do nothing (silent)
    elseif model.safety==2, model.error = @warning; % print warning
    elseif model.safety==3, model.error = @error;   % raise error
    end

    if model.safety>2
        % We will move to specifying linear reactions as a matrix
        % These are checked here for backwards compatibility / migration
        if isfield(model,'rQA'  ), warning('model.rQA is defined but not used'); end
        if isfield(model,'rAR'  ), warning('model.rAR is defined but not used'); end
        if isfield(model,'rRQ'  ), warning('model.rRQ is defined but not used'); end
        if model.sigma*model.n<0.4,
            model.error('Interaction radius finer than grid resolution');
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parse the model.description matrix
    if isfield(model,'linearRates') && isfield(model,'linearTransitions'),
        if model.verbosity>0 && isfield(model,'description'),
            'defined; using the latter.\n']);
        end
        model.linearRates  = model.linearRates(:);
        model.nstates      = size(model.linearTransitions,1);
        model.ntransitions = size(model.linearTransitions,2);
        assert (size(model.linearRates,1)==model.ntransitions);
    elseif isfield(model,'description'),
        if model.verbosity>0,
            if isfield(model,'linearRates'),
                'but missing linearTransitions. Using description matrix\n']);
            end
            if isfield(model,'linearRates'),
                'but missing linearRates. Using description matrix\n']);
            end
        end
        model.ntransitions      = size(model.description,1);
        model.nstates           = size(model.description,2)-1;
        model.linearRates       = model.description(1:end,model.nstates+1);
        model.linearTransitions = model.description(1:end,1:model.nstates)';

    else
        error('Please specify either model.description or model.linearRates and model.linearTransitions');
    end

    excitatory_reaction = [-1 1 zeros(1,model.nstates-2)]';
    reverse_excitatory_reaction = [1 -1 zeros(1,model.nstates-2)]';
         all(model.linearTransitions(1:end,1)==reverse_excitatory_reaction)),
            'reaction represents the Q to A transition.\n\n']);
        error('Invalid reaction description');
    end

    % Depends on model.linearRates
    model.Rlin         = model.linearTransitions*diag(model.linearRates)*(model.linearTransitions<0)';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % auxiliary variables
    model.dx     = model.n.^-2;      % Volume of each basis element (grid square)
    model.nn     = model.n*model.n; % Total number of spatial basis elements
    model.Inn    = eye(model.nn);   % Identity matrix for single-population
    model.volume = model.dx * abs(model.dt);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize observation model
    % Several variables interact to generate observed (expected) spikes

    % model.alpha is a dispersion parameter to account for non-poisson
    % firing. It is equal to one-over-the-square-root of the fano factor.
    % Poisson processes will have a fano factor of 1 (and an alpha of 1)
    % Rhythmic spiking will have a fano factor less than 1, and an alpha
    % greater than one. Irregular bursting can have a fano factor larger
    % than one, and an alpha less than one.
    % Latent rates should be multiplied by alpha (and the gain)

    if ~isfield(model,'inhomogeneous_gain'),
        model.inhomogeneous_gain = ones(model.n.^2,1);
    end
    if any(model.gain<0.0),
        error('All observation gain parameters must be positive');
    end
    if any(model.inhomogeneous_gain<0.0),
        error('All observation gain parameters must be positive');
    end

    model.goodchannels = (model.gain.*model.inhomogeneous_gain)>0.0;
    model.bias = model.bias(:) .* ones(model.n.^2,1);
    model.bias(~model.goodchannels) = 0.0;
    if any(model.bias<0.0),
        error('All observation bias parameters must be positive');
    end

    % Combine unform and spatially-varying gain parameters
    netgain = model.gain * model.inhomogeneous_gain(:);
    % Adjust the biases,
    adjbias = model.bias(:).^model.gamma;
    % and rescale by gains
    bias    = adjbias ./ (netgain * model.volume);
    % Multiply rates by spatiotemporal volume and dispersion parameter
    gain    = netgain * model.volume;% * model.alpha;
    bias(~isfinite(bias))=0.0;
    gain(~isfinite(gain))=0.0;
    model.adjusted_gain = gain;
    model.adjusted_bias = bias;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % precompute operator for nonlocal interactions
    minsigma = 1e-2/model.n;
    if model.sigma<minsigma,
        fprintf(2,'The interaction scale sigma is too small relative to the grid resolution\n');
        fprintf(2,'Implementing sigma %0.2e<%0.2e as a decoupled system\n',model.sigma,minsigma);
        if model.cutoff,
            warning('Cannot use minimum scale cutoff in a system with no interaction scale; disabling');
            model.cutoff==false;
        end
        % replace interaction kernels with dummy (identity) kernels
        model.K2D = eye(model.nn);
        model.K1D = eye(model.n);
    else
        model.K2D = gaussian2DblurOperator(model.n,model.sigma*model.n);
        model.K1D = gaussian1DblurOperator(model.n,model.sigma*model.n);
    end
    model.blur = @(x) reshape(model.K1D*reshape(x,model.n,model.n)*model.K1D',1,model.nn);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % These operators are used for faster updates in the moment closure
    % code. They extract the quiescent and active subspaces
    I  = model.Inn;
    O  = zeros(model.nn,model.nn);
    model.getQ = [I O];
    model.getA = [O I];
    for i=1:model.nstates-2,
        model.getQ = [model.getQ O];
        model.getA = [model.getA O];
    end
    model.getQ  = model.getQ';
    model.getA  = model.getA';
    model.getKA = (model.K2D*model.getA')';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Matrices defining important subspaces of the state space
    % Number of individual concentrations to model
    model.dimension = model.nstates*model.nn;
    % Identity matrix for full model
    model.diag = eye(model.nstates*model.nn);
    % Subspace reflecting spatial correlatoins and iteractions
    model.K2D3 = kron(eye(model.nstates),model.K2D);
    % Noise correlations for each reaction transitoin channel
    model.K2DT  = kron(eye(model.ntransitions),model.K2D);
    % Subspace orthogonal to local neuron count number and spatial correlations
    model.state_var_reg = model.K2D3*kron(eye(model.nstates)*model.nstates-ones(model.nstates),model.Inn)*model.K2D3';
    % Subspace orthogonal to spatial correlations
    model.count_var_reg = model.K2D3*model.K2D3';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Prepare matrices used for regularization
    % There are model.three subspaces that can be adjusted separately
    %
    % (1) Conserved regularization: `state_var_reg`
    %   Preserves spatial correlation and neuron density. Adds uncertainty
    %   only about which state neurons are in.
    %
    % (2) Non-conserved regularization: `count_var_reg`
    %   Preserves spatial correlations but not neuron density.
    %
    % (3) Diagonal regularization: `diag_var_reg`
    %   Reduces spatial correlations. Diagonal (shrinkage) regularization on
    %   covariance.
    %
    % Additive noise regularization matrix (and its Cholesky factor)
    % Adds variance only in directions orthogonal to the population count
    Reg = model.reg_state_var*model.state_var_reg;
    % Adds diagonal variance which incorporates uncertainty about the total count
    Reg = Reg + model.reg_count_var*model.count_var_reg;
    % Adds diagonal (purely independent) regularization
    Reg = Reg + model.reg_diag*model.diag;
    % Square-root update uses cholesky factor of regularization noise
    Reg = 0.5*(Reg+Reg');
    % Store pre-computed cholesky factor of regularization
    % (used in square-root form update)
    try
        model.regMsqrt = chol(Reg); % x=chol(x)'*chol(x)
    catch
        % Regulaization singular? Or zero?
        model.regMsqrt = Reg;
    end
    model.reg      = Reg;
    % Intermediate value used in covariance update
    model.Ireg     = model.diag + model.reg;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define model reaction system

    model.Sspatial = kron(model.linearTransitions,eye(model.nn));
    % Depends on model.linearRates and model.Rlin
    model.Rspatial = kron(model.Rlin,eye(model.nn));

    % model.noiseCov incorporates
    % 1) anti-correlations caused by reactions from one state to another
    % 2) spatial correlation from basis overlap (if using scale cutoff)
    % 3) scaling by spatial volume of basis functions
    % 4) scaling by dt "temporal volume" of integration step size
    % If using the cutoff option, the spatial basis functions are
    % overlapping gassians determined by the nonlocal interaction kernel
    % the full kernel for all states is stored in model.K2DT
    % The noise model assigned fluctuations from reactions to covariances
    % with correlations determined by the stoichiometry, and with
    % additional correlations induced by overlap in the basis functions.
    % We precompute this matrix product here.
    if model.cutoff,
        % Effective system size scalar if we set our minimum spatial scale to match
        % the interaction scale for the system. One advantage of using this
        % convention is that changing the simulation grid scale will not change
        % the underlying system being simulated. In this case the simulation grid
        % introduces only discretization error, but the underlying statistical
        % field theory is identical. If on the other hand we allow the simulation
        % grid to set the minimum scale for the model, then finer simulation grids
        % will add higher frequencies into the model. The statistical field theory
        % is not valid at scales below the size of single neurons, so too-fine a
        % a discretization may lead to invalid behavior. Using a fixed model scale
        % set by the excitatory interactions avoids this ambiguity.
        assert(model.sigma>minsigma);
        model.ss = 2*pi*model.sigma^2.*model.ss_rescale;
        model.noiseCov = model.Sspatial*model.K2DT;
    else,
        % Effective system size scalar if minimum scale is set by the simulation
        % grid. (This should typically be finer than the scale of interactions,
        % so simulating at this scale is like adding some more high-energy, i.e.
        % high spatial frequency, modes into the system. This increases the
        % apparent fluctuations. The statistical field theory becomes invalid
        % at spatial scales comparable to the size of single neurons, so using
        % the LNA at this finier spatial scale is expected to be less accurate.
        % This corresponds, effectively, to using a smaller population size,
        % which naturally has a larger error term in the system-size expansion)
            model.ss = model.dx.*model.ss_rescale;
            model.noiseCov = model.Sspatial;
            model.noiseCov = sparse(model.noiseCov);
    end
    % Adjust noise to reflect spatiotemporal basis function volume
    % We approximate the noise as the integral over space and time
    % of a Poisson random measure, with intensity proportional to the
    % reaction rates. `ss` is the system size in the expression below.
    % note: ss may need to be replaced by a vector if
    % basis functions of different volumes are used
    %
    % The absolute-value surronding `dt` is there to support (in the future)
    % negative time steps. Reversing time for the mean dynamics, while keeping
    % positive time for the fluctuations, in a simple way to approximate the
    % adjoint operator for the system, which allows for backward filtering.
    model.noiseCov = model.noiseCov.*sqrt(abs(model.dt)*model.ss);

    % Linear contribution to Jacobain (constant in time)
    % Depends on model.linearRates and model.Rlin
    model.Jlin  = kron(model.Rlin,eye(model.nn)).*model.dt;
    % An intermediate calculation in the covariance update
    model.IJlin = model.diag - model.Jlin;
    % Linear reactions can be integrated forward exactly via exponentiation
    % Depends on model.linearRates and model.Rlin
    model.Flin = expm(model.Rlin*model.dt);
    % Linear reactions can be integrated forward exactly via exponentiation: all spatial locations
    model.FlinSpatial = kron(model.Flin,model.Inn);
    model.FlinReactionSpatial = kron(model.linearTransitions*diag(model.linearRates),model.Inn)';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Construct initial covariance state
    C =     model.ini_state_var * model.state_var_reg;
    C = C + model.ini_count_var * model.count_var_reg;
    C = C + model.ini_reg_diag  * model.diag;
    C = 0.5*(C+C');
    try
        model.Pini = cinv(C);
    catch
        display('System is too large, increase the reg_diag parameter');
        error('matrix operations suffering from precision loss');
    end
    model.Cini = C;

    if ~isfield(model,'names'),
        for i=1:model.nstates,
            model.names{i} = sprintf('State%d',i);
        end
    else
        if size(model.names,2)~=model.nstates,
                size(model.names,2),model.nstates));
    end

    if ~model.cscale,
        model.cscale=ones(model.nn*model.nspecies,1);
    end

end




    '''#STOPCODE
