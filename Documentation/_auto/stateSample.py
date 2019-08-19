#!/usr/bin/env python
# -*- coding: UTF-8 -*-
def sampleSates(model,varargin):
    r'''

    The function `sampleSates(model,varargin)` will sample the evolution
    of the intensitiy fields and corresponding point-process observations
    from a given model.

    This does not support spatially inhomogeneous gain parameters.
    This does not support the dispersion parameter alpha. (uses alpha=1)

    This uses the complex Langevin equation (Schnoerr & al. 2014), so it
    may return negative intensities. Negative intensities are clipped to
    0 for the purpose of point-process sampling.

    Parameters
    ----------
    model : struct
        pre-initialized model struct; see `initialiezModel`

    Other Parameters
    ----------------
    Nsample : `int`, default 500
        Number of time-points to sample
    Nburn : `int`, default 200
        Number of initial burn-in samples (to approach steady-state)
    doplot : `bool`, default `false`
        Whether to plot states while filtering
    upscale : `int`, default 8
        Amount of upsampling to perform on shown arrays, if plotting
    skipsim : `int`, default 10
        Skip every `skipsim` ("skip-simulation") frames, if plotting
    save_figure : `bool`, default `false`
        Whether to save plots (if doplot is on) for animation
    effectivepopsize : `float`, default 50,
        Effective population size used
    oversample : `int`, default 10
        Spatiotemporal spiking resolution (relative to n); Spikes are
        sampled as a Poisson process using the active (A) intensity field
        upsampled by this factor.
    efraction : `float`, default 0.05
        Fraction of cells within local region that are excited by
        spontaneous wave events.

    Returns
    -------
    ini : `matrix`,
        Packed initial conditions for mean concentrations. Species are
        packed in the order defined by the model, where the second
        species reflects the active "A" state.
    xydata : `cell`, 1×Ntimes
        Cell array of spiking events. For each timepoint, we have a
        Nspikes×2 matrix, where column 1 is the x location and column 2
        is the y location of each spike, in the [0,1]² box.
    rates: `cell`, 1×Ntimes
        Basis-projected spiking events, with rates normalized by
        spatiotemporal volume.
    simulatedM: `cell`, 1×Ntimes
        Each cell array contains a MN² by 1 matrix of packed states, where
        M is the number of states (species), and N is the size of the N×N
        spatial basis, with N² basis elements total. Species are packed in
        order, with the spatial values packed in Matlab's default
        (column major) order.

    '''
    pass#SKIPME
    r'''#STARTCODE

    %|===== PARAM ====|=VAL=|= NOTE =|
        'efraction'   ,0.05,'Fraction Quiescent cells excited in spontaneous events'
    }, varargin);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up model for sampling (requires some special settings)
    % Matlab does pass-by-reference with copy-on-write
    % Changes to the model struct here remain local
    model = initializeModel(model);

    % For sampling the Langevin equation, noise is the same as that of the LNA
    % Square-root form makes sampling fluctuations simpler
    model.method   = 'LNA';
    model.sqrtform = true;

    % Check for legacy model definitions
    if isfield(model,'rQA'),
        r1 = model.linearRates(1);
        if r1~=model.rQA,
            msg = 'model.rQA=%0.2f but model.linearRates(1)=%0.2f';
            warning(sprintf(msg,model.rQA,r1));
        end
    end

    % Effective event rate
    model.event_rate = model.linearRates(1).*model.dt*model.dx;
    if (model.event_rate<1e-12),
        warning('Event rate <1e-12; is this intended?');
    end
    if (model.event_rate>0.1),
        display('Warning: spontaneous rate too nigh to model as Poisson noise');
        display('Increase resolution or reduce spontaneous rate?');
        model.event_rate = min(max(1.0-exp(-model.event_rate),0.0),1.0);
    end

    % For simluation, spontaneous depolarization is modeled as shot noise
    % which has a gain-like effect of increasing the Q->A transition in a
    % local area. The effect must be rescaled for neuron effectivepopsize when the
    % simulation grid scale is adjusted. This parameter not used for inference.

    % Old behavior
    % model.event_rescale = (model.n/20.0).^2;
    % To reproduce old behavior with new parameters, set
    % efraction = 0.04 ~= 0.25/(2*pi);

    % Wave events excite 5% of the quiescent cells within a local region
    model.event_rescale = 2*pi*(model.sigma.*model.n).^2.*opt.efraction;

    % Spontaneous excitation will be handeled separately.
    % Disable it here so it does not affect the mean update
    model.linearRates(1) = 0.0;

    % Re-initialize auxiliary variables to reflect new paramters
    model = initializeModel(model);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initial condition of intensity fields
    celleffectivepopsize = ones(1,model.n*model.n)*opt.effectivepopsize;
    Q = celleffectivepopsize*0.9;  % Quiescent
    A = celleffectivepopsize*0.05; % Active
    R = celleffectivepopsize*0.05/(model.nstates-2); % Refractory
    ini = [Q A kron(R,ones(1,(model.nstates-2)))];
    M = ini;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulate neural field model
    % Burn in to reach steady-state
    for i=1:opt.Nburn
        % Solve the model forward in time with known parameters
        M = meanUpdate(model,M);
        M = PoissonEvents(model,M);
    end

    % When plotting show all events since last displayed time point (so as not
    % to give the false impression that the algorithm can infer fields from
    % extremely sparse data)
    last_frame = 1;
    spatial_mean_projection = kron(eye(model.nstates),ones(1,model.nn))./model.nn;
    % Build a projection matrix for computing confidence intervals
    if opt.doplot,
        opt.figure = figure(17);
        NFCP_plotting.fix_figure(opt);
        NFCP_plotting.centerfig(gcf,1200,350);
    end

    % Remember the current state (initial conditions for sampling)
    ini = M;

    % Start recording (presumed at steady-state)
    xydata       = {}; % Store simulated spiking observations
    rates        = {}; % Store rate intensity field
    saved_images = {}; % Store neural field images for later display
    saved_means  = [];

    for i=1:opt.Nsample
        % sample spikes: interpolate field to get finer spatial resolution
        [xydata{i},rates{i}] = sampleSpikes(model,M,opt);
        simulatedM{i}        = M;

        % Store states and spatial means for later display
        means = spatial_mean_projection*M(:);
        saved_means(i,1:model.nstates+1) = [means; sum(means)];

        % Show plot of simulation state
        if opt.doplot && mod(i,opt.skipsim)==0
            NFCP_plotting.fix_figure(opt);
            newplot; hold off; subplot(131);
            trueRGB = NFCP_plotting.fieldsToRGB(abs(M),model.n,model.cscale,opt.upscale);
                sprintf('Simulating timepoint %d/%d',i,opt.Nsample));
            % Show point events (if any) since last displayed frame
            while last_frame<=i,
                if (size(xydata{last_frame},1)>0)
                    plot(xydata{last_frame}(:,1),xydata{last_frame}(:,2),'.w','MarkerSize',12);
                end
                last_frame = last_frame+1;
            end
            % Plot average intensities over time
            subplot(1,3,[2 3]);
            times = (1:size(saved_means,1))*abs(model.dt);
            plot(times,real(saved_means));
            ta = max(0,times(end)-100);
            xlim([ta,ta+100]);
            xlabel('Time (seconds)');
            ylabel('Mean concentration');
            legend([model.names 'Total']);
            hold off;
            if opt.save_figure,
                NFCP_plotting.save_figure(i,false,trueRGB,opt);
            end
            fr=getframe(gcf);
            clear fr;
            hold off;
        end


        % Solve the model forward in time with known parameters
        % Use linear noise approximation to generate fluctuations
        M = meanUpdate(model,M);
        L = noiseModel(model,M);
        noise = L*randn(model.nn*model.ntransitions,1);
        M = M + noise;
        M = PoissonEvents(model,M);
    end
end


    '''#STOPCODE
def PoissonEvents(model,M):
    r'''

    Add random depolarizatoins as shot noise (Poisson process)

    Parameters
    ----------
    model : struct
        pre-initialized model struct; see `initialiezModel`
    M : vector
        vector of mean intensities

    Returns
    -------
    M : vector
        updated vector of mean intentities with random excitation events
        applied

    '''
    pass#SKIPME
    r'''#STARTCODE
    N = model.nn;
    Q = M(1:N);
    A = M(N+1:2*N);
    blurred = reshape(model.blur(rand(model.n)<model.event_rate)*model.event_rescale,model.nn,1);
    events  = blurred(:).*Q(:);
    M(1:N)     = Q(:)-events;
    M(N+1:2*N) = A(:)+events;
end


    '''#STOPCODE
def sampleSpikes(model,M,opt):
    r'''

    Sample spikes from 2D intensity fields (uses the second field, by
    convention termed the `A` or active-state field).

    This does not support spatially inhomogeneous gain parameters.
    This does not support the dispersion parameter alpha.

    This actually does Bernoulli sampling, but on a spatially upsampled
    domain where the Poisson process can be approximated as Bernoilli.

    The purpose of this is to sample a single spike at each spatial
    location.

    Parameters
    ----------
    model : struct
        pre-initialized model struct; see `initialiezModel`
    M : vector
        vector of mean intensities

    Returns
    -------
    spikes : 2d array
        Array of spiking (x,y) locations

    '''
    pass#SKIPME
    '''#STARTCODE
    N = model.nn;
    A = M(N+1:2*N);
    if strcmp(model.link,'linear'),
        x = A;
    elseif strcmp(model.link,'log'),
        x = exp(A);
    elseif strcmp(model.link,'squared'),
        x = A.^2;
    end
    % volume  = model.dt.*model.dx.*model.gain;
    % rate    = volume.*max(0,real(A));
    gain    = model.inhomogeneous_gain.*(model.gain*model.volume);
    rate    = gain.*max(0,real(A));
    upfield = imresize(reshape(rate,model.n,model.n),opt.oversample);
    bias    = model.bias.*model.volume;
    if all(size(model.bias)==[1 1]),
        upfieldadj = upfield+bias
    else,
        bias       = imresize(reshape(bias,model.n,model.n),opt.oversample);
        upfieldadj = upfield+bias;
    end
    [y,x]  = find(rand(size(upfieldadj))<upfieldadj.*opt.oversample.^-2);
             (rand(size(x,1),2)-0.5)./(model.n*opt.oversample),1),0);
end



    '''#STOPCODE
