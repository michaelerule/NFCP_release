#!/usr/bin/env python
# -*- coding: UTF-8 -*-
def covarianceUpdate(model,M,P):
    r'''

    Covariance update for spatiotemporal system, operating on the
    precision (inverse covariance) matrix.

    This update is obsolete. The subspace measurement update now returns
    the covariance, instead of the precision.

    Parameters
    ----------
    P: positive definite matrix
        Full precision (inverse covariance) matrix for the
        system. The first two states should be the Quiescent (Q) and
        actively firing (A) states. Each species has n x n spatial grid
        packed in row-major order
    M: vector
        Means for concentrations of species Q, A, R (concatenated)
    model: struct
        Pre-initialized model struct. See `initializeModel`.

    Returns
    -------
    P : matrix
        Updated *precision* (inverse covariance) matrix.

    '''
    pass#SKIPME
    '''#STARTCODE

    warning('the covarianceUpdate function is deprecated');

    I = model.diag;
    J = jacobian(M,model);

    SIGMA = noiseModel(model,M,P)

    % Use Cholesky to compute covariance update without explicit inversion
    % This also operates on the square-root form for improved stability
    % The expression approximated by the commands below is covariance
    % P = inv(F inv(P) F' + Q)
    % F is the Jacobian F = exp(J)
    RB   = chol(P+I.*model.reg_precision)*(I-J);

    if model.sqrtform,
        % Square-root form assumes that regM is passed as Cholesky factor
        % To use same regularization as full form, use the following
        SRB = RB*SIGMA;
        RRB = RB*model.regMsqrt';
        SQ = chol(I+SRB*SRB'+RRB*RRB')'\RB;
    else,
        % Covariance noise, (not square-root form)
        SQ = chol(I+RB*(SIGMA+model.reg)*RB')'\RB;
    end
    P = SQ'*SQ;

    % Safety flag
    %   0 = Do not test at all
    %   1 = Test and silently patch unexpected numerical situations
    %   2 = Test and warn on unexpected numerical situations
    %   3 = Test and Hard-fail on unexpected numerical situations
    if model.safety>0,
        if model.safety<3,
            % Safety levels 1 and 2 will try to fix the problem
            repaired = false;
            while rcond(P)<eps,
                P = P + model.reg_precision.*model.diag;
                repaired= true;
            end
            % Safety level 2 will issue a warning that a problem was patched
            if repaired && model.safety==2,
                warning('Repaired ill-conditioned precision matrix');
            end
        end
        % Safety level 3 will hard-fail if there is a problem
        assertMatrixCondition(P,model.safety,0);
    end





    '''#STOPCODE
