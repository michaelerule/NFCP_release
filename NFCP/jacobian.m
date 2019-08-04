function J = jacobian(M,model)
    %{
    Returns the locally-linearized system for propagating the covariance
    term forward in time.

    Parameters
    ----------
    M: vector
        Means for concentrations of species Q, A, R (concatenated)
    model: struct
        Pre-initialized model struct. See `initializeModel`.
        
    Returns
    -------
    J : matrix
        System jacobian (linear and nonlinear), scaled by `model.dt`
    %}
    
    % Number of spatial basis functions
    N  = model.nn;
    % Coupling strengths for nonlocal <QA> interaction 
    K  = model.K2D;
    % Excito-excitatory rate parameter
    rE = model.rAA;
    % Time step
    dt = model.dt;
    % Jacobian contribution from nonlinear, nonlocal reactions
    % This assumes that the first state is Quiescent 
    % and the second state is Active
    Q   = M(1:N);
    A   = M(N+1:2*N);
    dQQ = - rE.*(diag(K*A(:)).*dt);
    dQA = -(rE.*(diag(Q).*dt))*K;
    % Jacobian contribution from linear reactions
    % This is constant in time and is pre-computed in the
    % initializeModel function
    J = model.Jlin;
    J(1:2*N,1:2*N) = J(1:2*N,1:2*N) + [dQQ  dQA; -dQQ -dQA];



