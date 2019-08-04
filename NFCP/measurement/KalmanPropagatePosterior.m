function [M2,C2,rPa,rCa,rMa] = KalmanPropagatePosterior(M,C,Pa,Ma,pPa,pMa)
    %{
    Propagate a measurement observation on the A state to the other latent
    states that are not directly observed.     
    
    This version uses an estimated posterior on the active state A to
    update the other states (Q,...). 
    
    For a version that uses a Gaussian measurement of A, see
    `KalmanPropagateObservation`
    
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
    
    Returns
    -------
    M2 : vector
        posterior means for full model (Q,A,...)
    C2 : positive definite matrxi
        posterior covariances for full model (Q,A,...)
    rPa : positive definite matrix
    rCa : positive definite matrix
    rMa : vector
    %}
    
    m   = length(Ma);
    
    % Generate a surrogate gaussian measurement to propagate update to
    % full joint model including non-observed states. This can be done
    % in one step since the non-observed states are independent of the
    % measurements condtitioned on the observed (active) state. 
    rPa   = pPa - Pa;          % Get information added by update, then
    chrPa = chol(rPa);         % Invert using cholesky factorization
    chrCa = (chrPa\eye(m))';   % to find the Gaussian covariance of 
    rCa   = chrCa'*chrCa;      % measurement update.
    rMa   = chrPa\(chrPa'\(pPa*pMa - Pa*Ma)); % Get measuremend mean. 
    
    Ca  = C(m+1:2*m,m+1:2*m);  % Active state slice of prior covariance
    Ma  = M(m+1:2*m,1);        % Active state slice of prior mean
    CCa = rCa + Ca;            % Kalman covariance update subexpression
    Cx  = C(m+1:2*m,1:end);    % subspace slice of prior covariance
    Kg  = (CCa\Cx)';           % Kalman gain
    M2  = M + Kg*(rMa - Ma);   % posterior mean
    C2  = C - Kg*Cx;           % posterior covariance
    M2  = M2(:);               % ensure column vector
    
    C1 = 0.5*(C2+C2');         % ensure symmetric
