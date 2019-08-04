#!/usr/bin/env python
# -*- coding: UTF-8 -*-
def filteringUpdateCovariance(model,M,C):
    r'''

    Combined mean and covaraince update for filtering using covariance
    matrix input / output (as opposed to prexision matrix).

    Parameters
    ----------
    model: struct
        Pre-initialized model struct. See `initializeModel`.
    M: vector
        Means for intensities of states. The first states must represent
        Quiescent (Q) neurons. The second states must represent actively
        firing (A) neurons. The remaining states can implement arbitrary
        linear dynamics (with Poisson noise).
    C: matrix
        Covariance matrix; Optional if using LNA, in which case covariance
        is not returned. (expects one output argument). If model.sqrtform
        is true, then C should be a cholesky factor of the covariance.

    Returns
    -------
    M : vector
        Updated state estiamte for mean intensities
    P : matrix
        Updated *precision* (inverse covariance) matrix.

    '''
    pass#SKIPME
    '''#STARTCODE

    N  = model.nn;
    I  = model.diag;
    K  = model.K2D;
    rE = model.rAA;
    dt = model.dt;
    m  = model.n.^2;

    % Experimental: turn of reactions at low rates?
    %MINR     = 1d-2;
    %aresmall = M<MINR;
    %small    = M(aresmall);
    %M(aresmall) = 0.0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MEAN UPDATE

    % Locally linearize the system and solve forward using exponentiation
    % This enforces neuron count conservation exactly. Local, linear
    % reactions are time homogenous and the discrete-time forward operator
    % can be computed in advance via exponentiation.
    % The locally-linearize nonlinear excito-exctatory interaction
    % is state-dependent and spatially inhomogenous.
    % Integrating using exponentiation is not feasible. However, if
    % the nonlinear reactions are integrated "separately", then there is
    % a quick closed-form for computing the conserved update.

    % The linear forward operator includes time-step dt and is integrated
    % exactly using exponentiation and saved with the model.
    M = model.FlinSpatial*M(:);

    % Nonlinear reations are handled as a second integration step
    %
    % Implementation note:
    %
    % The first reaction channel is hard-coded to be Q->A, and should
    % always be present. If this reaction is not needed, set the rate to
    % zero.
    Q  = M(1:N);
    A  = M(N+1:N*2);
    KA = K*A(:);
    EINTERACT = KA.*Q(:);

    % Moment closure calculates the expected <Q(K*A)> concentration
    % Using the system covaraince matrix. This can be computed by
    % considering projections of the system onto just the Q subspace
    % and onto the KA subspace. The code below uses the cholesky factor
    % of the precision matirx to simultaneously approximately invert
    % the precision matrix and perform these projections.
    if strcmp(model.method,'momentClosure'),
        % Compute correlation effects based on mean intensities
        % AFTER applying linear transition update. (The operators
        % should be applied as a product, i.e. one-after-another)
        %EINTERACT = EINTERACT + diag(model.getKA' * C * model.getQ);
        Cqa = C(m+1:2*m,1:m)';
        EINTERACT = EINTERACT + sum(model.K2D.*Cqa,2);
    end

    % Multiply the <QA> interaction by the exctiation rate
    ERATE = EINTERACT.*rE;

    % For updating the means, we include a fixed threshold. The nonlinear
    % reactions are locally linearized and converted to an effective rate
    % which is then integrated forward using exponentiation, to enforce
    % conservation
    rate  = ERATE./(Q(:)+model.minrate) - model.thr;
    rate(rate<0.0) = 0.0;
    delta = (exp(-rate*dt)-1).*Q(:);

    % Add excitatory changes to the new mean intensity
    Q          = Q(:)+delta;
    A          = A(:)-delta;
    M(1:N)     = Q(:);
    M(N+1:2*N) = A(:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COVARIANCE UPDATE

    % Precision-matrix based updates are useful when working with the
    % Laplace update, which operates on the precision matirx.
    % However, the Kalman-style update can be computed using covariances.
    % Here, we propagate the evolution of the covariance matirx
    % forward.

    % The covariance update requires the system Jacobian. The mean update
    % requires strict conservation, but the covariance update is more
    % tolerant: Euler integration can be used for speed. The contribution
    % of linear, local, homogeneous state transitions to the jacobian
    % is computed ahead of time and stored with the model as Jlin. Only
    % the nonlocal, nonlinear contributions need to be computed on every
    % time-step, and these are handled by local linearization.
    se = (dt*rE);

    % Recompute KA reflecting updated intensities
    KA = K*A(:);

    % Compute nonlinear contribution to jacobian
    dQQ = diag(se.*KA);
    dQA = bsxfun(@times,K',se.*Q')'; %=diag(se.*Q)*K;
    Jnonlin = zeros(model.dimension);
    Jnonlin(1:2*N,1:2*N) = [-dQQ -dQA; dQQ dQA];

    % See initializeModel for definition of FlinReactionSpatial
    RATES      = model.FlinReactionSpatial*M(:);
    RATES(1:N) = RATES(1:N) + ERATE;
    RATES      = max(model.minrate,RATES);
    noise = model.noiseCov*bsxfun(@times,model.noiseCov',RATES) + model.reg;

    F = I + model.Jlin + Jnonlin;
    assertPSD(C,model.safety);
    C = F*C*F';
    %assertPSD(C,model.safety);
    C = C+noise;
    assertPSD(C,model.safety);

    % Constrain model if normalized flag is set s.t. densities sum to 1.
    if model.normalized,
        M       = reshape(M(:),model.nn,model.nstates);
        total   = sum(M,2);
        rescale = 1./total;
        % Adjust means
        rescale = kron(rescale,ones(1,model.nstates));
        M = M.*rescale;
        M = M(:);
        % Adjust precision matrix to reflect rescaling
        rescale  = reshape(rescale,model.nn*model.nstates,1);
        rescaleC = diag(rescale);
        assertPSD(C,model.safety);
        C = rescaleC * C * rescaleC;
        assertPSD(C,model.safety);
    end

    C = 0.5*(C+C');
    assertPSD(C,model.safety);

    % Add back in reactions? from small rates?
    %M(aresmall) = M(aresmall) + small;





    '''#STOPCODE
