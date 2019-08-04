function [M2,C2,pPa,pMa] = KalmanPropagateObservation(M,C,Pa,Ma,rPa,rMa,rCa)
    %{
    Propagate a measurement observation on the A state to the other latent
    states that are not directly observed.     
    
    This version uses a Gaussian measurement of the active state A to
    update the other states (Q,...). 
    
    For a version that uses a Gaussian posterior estimate of A, see
    `KalmanPropagatePosterior`
    
    Parameters
    ----------
    M  : full joint prior mean packed as Q, A, ...
    C  : full joint prior covariance packed as Q, A, ...
    Ma : ptior mean for active A state
    Pa : prior precision (inverse covariance) for active A state
    pPa : matrix
        Posterior precision for A state, can be none (`[]`) if `rCa` and 
        `rMa` are provided.
    pMa : vector
        Posterior mean for A state, can be none (`[]`) if `rCa` and `rMa` 
        are provided.
    rCa : matrix
        Covariance matrix for measurement of A state. Can be none (`[]`)
        if both `pPa` and `pMa` are provided.
    rMa : matrix
        Mean vector for measurement of A state. Can be none (`[]`)
        if both `pPa` and `pMa` are provided.
    
    Returns
    -------
    M2 : vector
        posterior means for full model (Q,A,...)
    C2 : positive definite matrxi
        posterior covariances for full model (Q,A,...)
    pPa : positive definite matrix
        Posterior precision matrix (inverse covariance) for active state
    pMa : vector
        Posterior means for active state
    %}
    
    m   = length(Ma);
    
    % Compute subspace posterior from surrogate likelihood
    % (needed for estimating measurement likelihood)
    pPa   = rPa + Pa;
    chpPa = chol(pPa+RIa);
    pMa   = chpPa\(chpPa'\(rPa*rMa + Pa*Ma));
    
    Ca  = C(m+1:2*m,m+1:2*m); % Active state slice of prior covariance
    Ma  = M(m+1:2*m,1);       % Active state slice of prior mean
    CCa = rCa + Ca;           % Kalman covariance update subexpression
    Cx  = C(m+1:2*m,1:end);   % subspace slice of prior covariance
    Kg  = (CCa\Cx)';          % Kalman gain
    M2  = M + Kg*(rMa - Ma);  % posterior mean
    C2  = C - Kg*Cx;          % posterior covariance
    M2  = M2(:);              % ensure column vector
