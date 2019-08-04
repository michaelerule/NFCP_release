function QA = excitatoryRate(model,M,P)
    %{
    Computes the expected product of the intensities for the quiescent and
    active states. Used to estiamte the rate of excitation. 
    
    If moment closure is enabled, the covariance is used to get a more
    accurate estimate. This requires inverting the precision matrix. 
    If the precision is close to singular, moment closure may fail. The
    parameter `model.reg_inverse` controls how much to regularize when
    inverting to recover the covariance. 
    
    Moment closure for more general functions (e.g.) a firing-rate 
    nonlinearity may be approximated by a locally quadratic function
    (if fluctations remain small). 
    
    Parameters
    ----------
    model: struct
        Pre-initialized model struct. See `initializeModel`.
    M: vector
        Means for intensities of states. The first states must represent
        Quiescent (Q) neurons. The second states must represent actively
        firing (A) neurons. The remaining states can implement arbitrary
        linear dynamics (with Poisson noise). 
    P: matrix
        Precision matrix
        
    Returns
    -------
    QA : vector
        Expected product of Q and A (local)
    %}
    N  = model.nn;
    K  = model.K2D;
    Q  = M(1:N);
    A  = M(N+1:N*2);
    QA = (K*A(:)).*Q(:);
    if strcmp(model.method,'momentClosure'),
        
        % Regularize for inversion
        P  = P+model.diag*model.reg_inverse;
        % Invert using Cholesky
        C  = cinv(P);
        % Extract Q-A block and apply spatial interactions
        QA = QA + diag(K*C(N+1:N*2,1:N));
        
        %{
        % Here is a faster way
        ch    = chol(P+model.diag*model.reg_inverse);
        CKAQ5 = (K*model.getA')*(ch\(ch'\model.getQ));
        QA = QA + diag(CKAQ5);
        %}
        
        % No negative rates!
        QA(QA<0) = 0;
    end
    % Scale by excitatory rate parameter
    QA = QA.*model.rAA;

