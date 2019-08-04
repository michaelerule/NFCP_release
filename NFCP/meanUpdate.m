function M = meanUpdateM(model,M,P)
    %{
    Conserved state update for mean intensities
    
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
        Precision matrix; Optional if using LNA
    
    Returns
    -------
    M : vector
        Updated state estiamte for mean intensities
    %}
    if nargin<3,
        % No precision matrix provided
        if strcmp(model.method,'momentClosure'),
            error('Moment-closure update requires the precision matrix');
        end
        P = eye(model.dimension);
    end
    % Propagate linear transitions with exponentiation
    M = model.FlinSpatial*M(:);
    % Propagate nonlinear transitions            
    N = model.nn;
    Q = M(1:N);
    A = M(1*N+1:2*N);
    % Locally approximate nonlinear interactions as linear
    rate = excitatoryRate(model,M,P)./(Q(:)+model.minrate) - model.thr;
    % Ensure no negative rates
    rate(rate<0.0) = 0.0;
    % Integrate linear rates forward with conservation
    % Considering only this (linearized) reactoin, exponentiation solves
    % forwared while inforcing conservation and positivity
    delta = (exp(-rate*model.dt)-1).*Q(:);
    % Update intensities
    M(1:N)       = Q(:)+delta;
    M(1*N+1:2*N) = A(:)-delta;
