#!/usr/bin/env python
# -*- coding: UTF-8 -*-
def filteringUpdate(model,M,P):
    r'''

    Combined mean and covaraince update for filtering.
    This operates on the precision matrix.

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
    P : matrix
        Updated *precision* (inverse covariance) matrix.

    '''
    pass#SKIPME
    '''#STARTCODE

    assertMatrixCondition(P,model.safety);
    N  = model.nn;
    I  = model.diag;
    K  = model.K2D;
    rE = model.rAA;
    dt = model.dt;

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
    %
    % The linear forward operator includes time-step dt and is integrated
    % exactly using exponentiation and saved with the model.
    M  = model.FlinSpatial*M(:);

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

    % Moment closure and noise updates work on covariance, while
    % measurement updates work on inverse covariance (precision)
    % For improves stability for those matrix inverses, we compute a
    % (regularized) cholesky factor of the precision matrix
    chP = chol(P+I.*model.reg_precision);

    % Moment closure calculates the expected <Q(K*A)> concentration
    % Using the system covaraince matirx. This can be computed by
    % considering projections of the system onto just the Q subspace
    % and onto the KA subspace. The code below uses the cholesky factor
    % of the precision matirx to simultaneously approximately invert
    % the precision matrix and perform these projections.
    if strcmp(model.method,'momentClosure'),
        % Compute correlation effects based on mean intensities
        % AFTER applying linear transition update. (The operators
        % should be applied as a product, i.e. one-after-another)
        chPKA = chP'\model.getKA;
        chPQ  = chP'\model.getQ;
        % The expression sum((chPKA.*chPQ),1) is a faster way of
        % computing the diagonal diag(chPKA*chPQ)
        EINTERACT = EINTERACT + sum((chPKA.*chPQ),1)';
    end

    % Multiply the <QA> interaction by the exctiation rate
    ERATE = EINTERACT.*rE;

    % For updating the means, we include a fixed threshold. The nonlinear
    % reactions are locally linearized and converted to an effective rate
    % which is then integrated forward using exponentiation, to enforce
    % conservation of neuron numbers
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

    dQQ = diag(se.*KA);
    dQA = bsxfun(@times,K',se.*Q')'; %=diag(se.*Q)*K;
    Jnonlin = zeros(model.dimension);
    Jnonlin(1:2*N,1:2*N) = [-dQQ -dQA; dQQ dQA];

    % The covariance update operates on the covariance, but the measurement
    % update operates on the inverse covaraince (precision) matrix. The
    % following lines take advantage of various numerical tricks to update
    % the precision matrix without explicitly inverting it.

    % See initializeModel:  model.IJlin = model.diag - model.Jlin;
    RB = chP*(model.IJlin-Jnonlin);
    RATES = model.FlinReactionSpatial*M(:);

    % For efficiency, ERATE is taken from the concentrations after the
    % linear update and prior to the excitatory update. This may lead
    % to some errors, TODO check this.
    RATES(1:N) = RATES(1:N) + ERATE;

    % Scaling of the reaction noise RATES by the time step dt has been
    % absorbed into the pre-computed factor `noiseCov`
    S  = bsxfun(@times,model.noiseCov,sqrt(abs(RATES))');

    % If excitation is too large, this expression will generate large
    % eigenvalues making the precision matrix ill conditioned. Either
    % S has too-small eigenvalues, or RB has too-large eigenvalues.
    % S is related to the noise terms and, due to conservation, will be
    % singular. RB is the culprit. It has a mix of large and small
    % eigenvalues. The nonlinear contributions to the Jacobian affect
    % only a subspace, and when they have a very different magnitude
    % from the linear contributions, it causes RB to become ill-
    % conditioned.
    SQ = chol(I+RB*(S*S'+model.reg)*RB')'\RB;
    if model.safety>0,
        if rcond(SQ)<sqrt(eps),
            SQ = SQ + model.reg_precision.*model.diag;
            if model.safety>1,
                model.error('Precision matrix factor ill conditioned');
            end
        end
        assertMatrixCondition(SQ,model.safety);
    end
    P = SQ'*SQ;

    % Constrain model if normalized flag is set such that densities
    % sum to 1.
    if model.normalized,
        M       = reshape(M(:),model.nn,model.nstates);
        total   = sum(M,2);
        rescale = 1./total;
        % Adjust means
        rescale = kron(rescale,ones(1,model.nstates));
        M = M.*rescale;
        M = M(:);
        % Adjust precision matrix to reflect rescaling
        rescale = reshape(rescale,model.nn*model.nstates,1);
        rescaleP = diag(1./rescale);
        P = rescaleP * P * rescaleP;
    end

    assertMatrixCondition(P,model.safety);








    '''#STOPCODE
