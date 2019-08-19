function [loglikelihood,M2,C2,rMa,rPa,rCa] = measurementSubspace(model,M,C,xypoints,srMa,srPa,srCa)
    %{
    Laplace approximation of measurement update.
    
    This update is selected by setting `model.update` to 'Laplace-subspace'.
    
    This performs the update only on the subspace that interacts directly
    with the measurements. The rest is conditionally gaussian. Latent 
    activations that are coupled to the spiking observations are constrained
    to be positive. The rest of the states are treated as conditionally
    Gaussian, and may be set to negative means after propagating the 
    conditional update.
    
    If "surrogate" observations for the A state, provided as `srMa`,`srPa`,
    and `srCa`, reflecting the measured means, precisions (inverse
    covariance), and covariance are provided, then this function ignores
    the provided data `y` and intead performs the Kalman filter update
    using these Gaussian measurements. 
    
    The (optionally) returned values `rMa`,`rPa`, and `rCa` reflect the
    means, precisions (inverse covariance), and covariance (respectively)
    of a Gaussian approximation of effective measurement, which can be used 
    with the Kalman filter measurement update to recover the 
    same posterior, given the provided prior (M,C). Can be stored
    and re-used as a faster approximate measurement for small 
    changes in parameters.
        
    Parameters
    ----------
    model: struct
        Pre-initialized model struct. See `initializeModel`.
    M: vector
        Means for concentrations of species Q, A, R (concatenated)
    C: positive definite matrix
        Full covariance matrix for the
        system. The first two states should be the Quiescent (Q) and
        actively firing (A) states. Each species has n x n spatial grid 
        packed in row-major order
    xypoints : npoints x 2 matrix
        Point xypoints in the form of a list of (x,y) locations.
        All points should like in the [0,1]² unit square.
        The function will also behave correctly if xypoints have been
        binned beforehand to an n x n grid.
    srMa,srPa,srCa : 
        Optional: pass Gaussian observation means, precisions, and 
        covariances. If these are provided, `y` will be ignored.
        
    Returns
    -------
    loglikelihood : 
        The total log-likelihood of the xypoints after measurement 
        update, excluding constant factors, computed via Laplace
        approximation.
    M2 : vector
        the posterior mode concentration for each location
    P2 : positive definite matrix
        the posterior precision for each location
    rMa : vector
    rPa : positive definite matrix
    rCa : positive definite matrix
    %}

    if nargin==4, use_surrogates = false; end
    if nargin>4,
        if nargin~=7,
            error('If using surrogate likelihoods, mean, precision, and covariance should be provided')
        end
        use_surrogates = true;
    end
    
    assertPSD(C,model.safety);       % ensure is PSD
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract active state
    m    = model.nn;                  % No. basis functions per state
    M    = M(:);                      % ensure means as column vector
    Ma   = M(m+1:2*m,1);              % extract active state means
    Ca   = C(m+1:2*m,m+1:2*m);        % extract active state covariance
    assertPSD(Ca,model.safety);       % ensure active covariance is PSD
    RIa  = model.reg_inverse.*eye(m); % Define diagonal regularizer
    chCa = chol(Ca+RIa);              % Invert using cholesky decomposition
    chPa = (chCa\model.Inn)';         % x = chol(x)'*chol(x)
    Pa   = chPa'*chPa + RIa;          % active state inverse covariances
    assertPSD(Pa,model.safety);       % ensure is PSD
    
    if use_surrogates,
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Surrogate update: ignore y, use provided Gaussian observations
        rMa = srMa;
        rPa = srPa;
        rCa = srCa;
        % Kalman update to get joint from updated A state
        [M2,C2,pPa,pMa] = KalmanPropagateObservation(M,C,Pa,Ma,rPa,rMa,rCa);
    else,    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Subspace Laplace approximation update
        y = binCounts(model,xypoints);  
        [pMa,pPa] = laplaceUpdateNewtonRaphson(model,Ma,Pa,Ca,chPa,y);
        % Kalman update to get joint from updated A state
        [M2,C2,rPa,rCa,rMa] = KalmanPropagatePosterior(M,C,Pa,Ma,pPa+RIa,pMa);
    end 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimate of the likelihood
    if model.dolikelihood,
        loglikelihood = logLikelihood(model,M,cinv(C),M2,cinv(C2),y);
        assert(~isnan(loglikelihood));
        assert(all(isfinite(loglikelihood)));
    else
        loglikelihood = NaN;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optional negative-means handling (slow!)
    if any(M2<0),
        % Could re-run constrained Laplace update? 
        % [loglikelihood,pmode,pprec] = measurementLaplace(model,M,P,xypoints)
        % [loglikelihood,M2,C2,rMa,rPa,rCa] = measurementSubspace(model,M,C,xypoints,srMa,srPa,srCa)
        if model.verbosity>=2,
            fprintf(2,'Subspace update yielded negative intensities; using constrained Laplace.\n');
        end
        [loglikelihood,M2,P2] = measurementLaplace(model,M,cinv(C),xypoints);
        C2 = cinv(P2);
        %{
        % Global correction (considers joint covariance)
        % M2 = nearestPositivePoint(M2,cinv(C2+eye(model.dimension)*model.reg_inverse));
        % Local correction
        for location=1:model.nn,
            idxs = (0:model.nstates-1).*model.nn+location;
            M_local = M2(idxs);
            if any(M_local<0),
                C_local = C2(idxs,idxs);
                M_local = nearestPositivePoint(M_local,cinv(C_local+eye(model.nstates)*model.reg_inverse));
            end
            M2(idxs) = M_local;
        end
        %}
    end
    
    
    

        
    
    








