#!/usr/bin/env python
# -*- coding: UTF-8 -*-
def stateInfer(ini,model,xydata,true_states,varargin):
    r'''

    Infer latent states given point-process observations

    Note: "Surrogate" measurements: In the filtering algorithm, we use a
    multivariate Gaussian approximation for latent states. We estimate
    an approximate Gaussian posterior for the non-conjugate Poisson
    measurement updates. Estimating this approximation can be time
    consuming. Once an approximation is estimated, it is possible to
    divide out the Gaussian prior to generate a conjugate Gaussian
    "measurement" that has the same effect as the Poisson measurement.
    This Gaussian measurement can be re-used for nearby parameterizations,
    and is much faster than recomputing the approximation. We call these
    Gaussian measurements "surrogate" measurements.

    Parameters
    ----------
    ini : `array`
        Initial mean intensities
    model : `struct`
        pre-initialized model structure; see `initializeModel`
    true_states : `cell`, default `false`
        Ground-truth states (if avaialable).
        If unavailable, defaults to `false`.

    Other Parameters
    ----------------
    doplot : `bool`, default `false`
        Whether to plot states while filtering
    save_figure : `bool`, default `false`
        Whether to save plots (if doplot is on) for animation
    upscale : `int`, default 8
        Amount of upsampling to perform on shown arrays, if plotting is on.
    skipinf : `int`, default 10
        Skip every `skipinf` ("skip-inferred") frames if plotting is turned
        on.
    showduration : `int`, default 200
        Total number of time-steps to plot (if plotting)
    showmaxy : `int`, default 2
        Y-axis maximum for plotting inferred states
    showprog : `int`, default false
        Whether to print debugging information for every frame
    normll : `bool`, default true
        Whether to normalize log-likelihood by the number of data
        samples.
    timeout : `float`, default -1
        If this is a positive number, monitor seconds elapsed and terminate
        the algorithm early if the runtime exceeds `timeout` seconds.
        Note: not supported yet
    ratescale : `float`, default 1
        The latent A state will be multiplied by this value for display
    peakactivity : `bool`, default false
        If true, plot the active state density at the location with
        the most activity, rather than the spatial average
    mergeRstates : `bool`, default true
        Merge all refractory states when plotting
    poolpoints : `int`, default 5
        Number of frames of point events to plot. Should not exceed
        the displayed frame rate specified by `skipinf`
    MarkerSize : `int`, default 3
        MarkerSize parameter forwarded to scatter for plotting point
        events
    softmax : bool, default true
        If true, normalize the color-scheme to the peak value in each
        channels in plots.
    opt.figure : integer
        Restore focus to this window before each plotting call.
        this is a workaround for a Matlab bug, wherein a user focus event
        can cause new plots to be drawn to the wrong figure.
    iniP : matrix
        Optional initial conditions for correlation structure
    rescale : bool
        Whether to normalize plot output to reflect the fraction of cells in
        each state, as opposed to the absolute number of cells.

    Returns
    -------
    llsum : `float` or `vector`
        Log-likelihood of data given model (up to a constant factor),
        summed over all timepoints. Depending on model configuration,
        the log-likelihood may be decomposed into several contributions,
        in which case a vector is returned.
    infstate : `cell`, 1×Ntimes
        Cell array of mean intensity vectors over time. Each cell array
        contains a MN² by 1 matrix of pacted states, where M is the number
        of states (species), and N is the size of the N×N spatial basis,
        with N² basis elements total. Species are packed in order, with
        the spatial values packed in Matlab's default (column major) order.
    margvar : `cell`, 1×Ntimes
        Cell array of marginal variances for each state, location. This
        contains the marginal variances (for each spatial location
        and species). It is packed and organized similarly to the `infstate`
        variable that contains the mean concentrations. It corresponds
        to the diagonal of the covariance matrix.
    infe : `cell`, 1×Ntimes
        Cell array of effective nonlocal excitation over time. It contains
        the local mean-excitation value at all timepoints. Each timepoint
        contains a N²x1 matrix. This reflects the total excitatory input
        into each spatial region at each timepoint, multiplied by the number
        or density of quiescent (Q) agents in that region.

    '''
    pass#SKIPME
    r'''#STARTCODE
    %|===== PARAM ====|= VAL =|
        'doplot'       ,false, 'whether to show plot during inference'
        'save_figure'  ,false, 'whether to save figures to disk for animation'
        'cache'        ,false, 'whether to cache the results'
        'upscale'      ,8    , 'upscaling factor for plotting fields'
        'skipinf'      ,10   , 'number of timesteps to skip when plotting'
        'showduration' ,200  , 'time axis length; this parameter is obsolete'
        'showmaxy'     ,2    , 'y axis height'
        'showprog'     ,false, 'print progress in terminal'
        'normll'       ,true , 'normalize the log-likelihood by No. samples'
        'timeout'      ,-1   , 'whether to give up after a certain time'
        'ratescale'    ,1    , 'rescale rates'
        'peakactivity' ,false, 'normalize activity to peak state'
        'mergeRstates' ,true , 'merge refractory states when plotting'
        'poolpoints'   ,5    , 'show points from multiple frames'
        'MarkerSize'   ,3    , 'marker size for plotting point events'
        'points'       ,true , 'whether to draw point-process events'
        'softmax'      ,false, 'option for rescaling field colors'
        'usesurrogate' ,false, 'if true, pass surrogate likelihoods to Laplace-subpace'
        'getsurrogate' ,false, 'if true, get surrogate likelihoods from Laplace-subpace'
        'figure'       ,17   , 'figure number or handle to use; defaults to fig 17'
        'rescale'      ,true , 'normalize plot to reflect fractions rather than absolut population size'
        }, varargin);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Argument check
    [model,opt] = verifyArguments(ini,model,xydata,true_states,opt);
    Nsample = model.Nsample;
    if opt.usesurrogate,
        sM = opt.sM;
        sP = opt.sP;
        sC = opt.sC;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Use forward filtering to recover latent states and likelihood

    % Accumulate log-likelihood
    llsum = 0;

    % Initial conditions for means, and covariance/precision
    M = ini(:);
    P = opt.iniP;
    C = cinv(P);

    % Build a projection matrix for computing confidence intervals
    % This projects the covariance onto the marginals for each state
    % combining all spatial locations.
    % If opt.mergeRstates is set, then all refractory states (all states
    % after the first two) are merged into a single aggregate state.
    if opt.mergeRstates,
        varscales = [zeros(2,model.nstates); ones(1,model.nstates)];
        varscales(1:3,1:3) = eye(3);
    else,
        varscales = eye(model.nstates);
    end
    % The active state (state 2) can be rescaled, since it is typically
    % small, for plotting
    varscales(2,2)= opt.ratescale;
    spatial_means = kron(varscales,ones(1,model.nn));
    total_mean    = ones(1,model.dimension);
    opt.meanproj  = [spatial_means; total_mean]./model.nn;
    opt.K         = size(opt.meanproj,1);
    inf_means     = zeros(Nsample,opt.K+1);
    inf_errs      = zeros(Nsample,opt.K);
    saved_means   = zeros(Nsample,opt.K);

    % If there is ground-truth data available, generate
    % variables used for plotting inferred vs. ground tuth
    if opt.have_groundtruth,
        for i=1:Nsample,
            trueM = real(true_states{i});
            saved_means(i,1:opt.K) = opt.meanproj*trueM(:);
        end
    end

    % If plotting, show all events since last displayed time point
    if opt.doplot,
        opt.figure = figure(17);
        NFCP_plotting.fix_figure(opt);
        NFCP_plotting.centerfig(gcf,1200,350);
    end

    % Monitor elapsed time: report progress, erasing previous message
    ttime=clock; stime=ttime;
    if opt.showprog, fprintf('\nProcessing frame '); todel=0; end

    for isamp=1:Nsample
        % Monitor progress: print message indicating progress every .5 s
        if opt.showprog && etime(clock,ttime)>1.0,
            todel = overprint(sprintf('%d of %d',isamp,Nsample),todel);
            ttime = clock;
        end
        if opt.doplot,
            if ~opt.use_covariance,
                % Store spatial marginal variances
                C = cinv(P+model.diag.*model.reg_precision);
            end
            % Plot progress
                xydata,true_states,saved_means,inf_means,inf_errs,opt);
        end
        if opt.timeout>0 && etime(clock,stime)>opt.timeout,
            % Algorithm has exceeded its time limit!
            warning(sprintf('Algorithm timed out at frame %d;',isamp));
            warning('Likelihood is not valid');
            llsum = llsum*Nsample/(i-1);
            break;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Iterate filtering
        % Remember means before state update
        M0 = M(:);
        if opt.use_covariance,
            % Measurement update using covariance matrix
            [M,C,EINTERACT] = filteringUpdateCovariance(model,M,C);
            assertPSD(C,model.safety);
            % Check that rates are non-increasing if no spikes
            % observed
            % (ad-hoc stability measure)
            % [~,Mzero,~] = opt.measurementUpdate(model,M,C,zeros(0,2));
            % if ~all(Mzero(:)<=M0(:)),
            %     warning('Combined state update and measuement predicts increasing rates even without spikes! Consider reducing excitation.');
            % end
            if opt.getsurrogate,
                [ll,M,C,rM,rP,rC] = opt.measurementUpdate(model,M,C,xydata{isamp});
                sM{isamp} = rM;
                sP{isamp} = rP;
                sC{isamp} = rC;
            elseif opt.usesurrogate,
                    sM{isamp},sP{isamp},sC{isamp});
            else
                [ll,M,C] = opt.measurementUpdate(model,M,C,xydata{isamp});
            end
            assertPSD(C,model.safety);
        else
            % Measurement update using precision matrix
            [M,P,EINTERACT] = filteringUpdate(model,M,P);
            P = assertPSD(P,model.safety);
            % Check that rates are non-increasing if no spikes
            % observed
            % (ad-hoc stability measure)
            %[~,Mzero,~] = opt.measurementUpdate(model,M,cinv(C),zeros(0,2));
            %if ~all(Mzero(:)<=M0(:)),
            %    warning('Combined state update and measuement predicts increasing rates even without spikes! Consider reducing excitation.');
            %end
            if opt.getsurrogate,
                [ll,M,P,rM,rP,rC] = opt.measurementUpdate(model,M,P,xydata{isamp});
                sM{isamp} = rM;
                sP{isamp} = rP;
                sC{isamp} = rC;
            elseif opt.usesurrogate,
                    sM{isamp},sP{isamp},sC{isamp});
            else
                [ll,M,P] = opt.measurementUpdate(model,M,P,xydata{isamp});
            end
            assertPSD(P,model.safety);
        end
        llsum = llsum+ll;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Save filtered states;
        infstate{isamp} = M;
        infe{isamp} = EINTERACT;

        % Get the marginal variances
        if opt.use_covariance, margvar{isamp}=diag(C);
        % Otherwise we need to invert the precision matrix
        else,
            try, margvar{isamp}=diag(cinv(P));
            catch exception,
                % safety 0=unchecked, 1=patch, 2=patch and warn, 3=abort
                if model.safety==0, rethrow(exception);
                else,
                    P = P+model.reg_precision.*model.diag;
                    model.error(sprintf('Singular precision in frame %d',isamp));
                end
            end
        end
    end
    if opt.showprog, fprintf('\n'); end

    % Normalize l.l. by samples (in case timeout caused early termination)
    if opt.normll, llsum = llsum/isamp; end
    if model.verbosity>0,
        display('Finished filtering');
        display(sprintf('Total negative log-likelihood (un-normalized) is %0.4e',llsum));
    end
end


    '''#STOPCODE
def verifyArguments(ini,model,xydata,true_states,opt):
    r'''

    Argument verification

    '''
    pass#SKIPME
    r'''#STARTCODE
    model = initializeModel(model);
    % Sanity check number of samples given for inference
    Nsample = size(xydata,2);
    if isfield(model,'Nsample') && Nsample~=model.Nsample,
        error(sprintf('model.Nsample=%d, but %d samples given for inference?',model.Nsample,Nsample));
    else,
        model.Nsample = Nsample;
    end
    % Detect whether ground-truth states were provided
    opt.have_groundtruth = ~islogical(true_states)||true_states~=false;
    if (opt.usesurrogate|opt.getsurrogate),
        if ~strcmp(model.update,'Laplace-subspace'),
            'Surrogate measurements only supported with Laplace-subspace update');
        end
        if opt.usesurrogate==opt.getsurrogate,
            'Can''t use and compute surrogate likelihoods simultaneously');
        end
    end
    if nargout>4,
        if ~opt.getsurrogate,
            'To get surrogate likelihoods, use Laplace-subspace and getsurrogate=true');
        end
        if nargout~=7,
            'method can optionally return surrogate means, precisions, covariances');
        end
    end
    if nargout==7 & ~opt.getsurrogate,
        'Surrogate means, precisions, covariances are returned as three extra outputs');
    end
    if opt.usesurrogate,
        if ~isfield(opt,'sM'),
            'To use surrogate likelihoods, pass surrogate means as ''sM'',values');
        end
        if ~isfield(opt,'sP'),
            'To use surrogate likelihoods, pass surrogate precisions as ''sP'',values');
        end
        if ~isfield(opt,'sC'),
            'To use surrogate likelihoods, pass surrogate covariance as ''sC'',values');
        end
        sM = opt.sM;
        sP = opt.sP;
        sC = opt.sC;
    end
    % Initial conditions for means
    M = ini(:);
    if ~(numel(M)==model.nstates*model.n*model.n),
        'Initial state must be the same size state-space: Nspecies x Nspatial^2');
    end
    % Initial conditions for the precision matrix
    % These are set heuristically and stored with the model is initialized
    if ~isfield(opt,'iniP'),
        opt.iniP = model.Pini;
    end
    assertPSD(opt.iniP,model.safety);
    assertPSD(cinv(opt.iniP),model.safety);

    % Bind measurement update function based on parameters
    opt = getUpdateModel(model,opt);
end


    '''#STOPCODE
def getUpdateModel(model,opt):
    r'''

    Depending on which update is used,
    keep track of either the precision or covariance matrix

    Laplace update needs the precision matrix.

    Kalman update can use the covariance matrix directly, which has
    a faster update.

    '''
    pass#SKIPME
    '''#STARTCODE
    if strcmp(model.update,'Laplace'),
        opt.use_covariance = false;
        opt.measurementUpdate = @measurementLaplace;
    elseif strcmp(model.update,'Laplace-subspace'),
        opt.use_covariance = true;
        opt.measurementUpdate = @measurementSubspace;
    elseif strcmp(model.update,'Kalman'),
        error('The "Kalman" update is deprecated.');
        opt.use_covariance = true;
        opt.measurementUpdate = @measurementKalmanCovariance;
    elseif strcmp(model.update,'Gamma-subspace'),
        error('The Gamma moment-matching update is not implemented yet');
        opt.use_covariance = true;
        opt.measurementUpdate = @measurementGamma;
    elseif strcmp(model.update,'EP'),
        error('Expectation propagation is not supported yet.');
        opt.use_covariance = true;
        opt.measurementUpdate = @measurementEP;
    else,
        error('Update method must Laplace, Kalman, Laplace-subspace, or Gamma-subspace');
    end
end



    '''#STOPCODE
