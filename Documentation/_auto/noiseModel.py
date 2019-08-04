#!/usr/bin/env python
# -*- coding: UTF-8 -*-
def noiseModel(model,M,P):
    r'''

    Computes the covariance update for fluctuations based on Gaussian
    moment closure. This is combined with propagation of the
    covariance from the previous time steps, as well as regularizing
    noise, to produce the covariance state update for filtering.

    Parameters
    ----------
    model: struct
        Pre-initialized model struct. See `initializeModel`.
    M: vector
        Means for intensities of states. The first states must represent
        Quiescent (Q) neurons. The second states must represent actively
        firing (A) neurons. The remaining states can implement arbitrary
        linear dynamics (with Poisson noise).
    P: `matrix`
        Precision matrix.
        Only required if using `model.method='momentClosure'`

    Returns
    -------
    noise : `matrix`
        If `model.sqrtform=true`, this will be a square-root of the noise
        covariance matrix. This is useful for sampling: the product of
        this factor with a vector of i.i.d. Gaussian variables will generate
        noise with the appropriate covariance for sampling the Langevin
        equation. The state-inference code can also compute some
        expressions using this square-root for improved numerical accuracy.
        If `mode.sqrtform=false`, this will return the full covariance
        matrix. This is used in the state-inference filtering code.

    '''
    pass#SKIPME
    '''#STARTCODE
    if nargin<3,
        % No precision matrix provided
        if strcmp(model.method,'momentClosure'),
            error('Moment-closure update requires the precision matrix');
        end
        P = eye(model.dimension);
    end

    % Get number of spatial basis functions
    N = model.nn;
    % Linear, local rates are spatially homogeneous
    RATES = kron(model.linearTransitions*diag(model.linearRates),model.Inn)'*M(:);
    % Nonlinear, nonlocal rates
    RATES(1:N) = RATES(1:N) + excitatoryRate(model,M,P);
    % Integrate one time-step (independent, linear approximation)
    %RATES = RATES;%;.*model.dt;
    % Integration subsumed into model.noiseCov

    if ~(isfield(model,'complexlangevin') && model.complexlangevin),
        RATES(RATES<0) = 0;
    end

    if model.sqrtform,
        % Return the square-root of the noise covariance matrix
        % This is useful for sampling and for square-root filters
        noise = model.noiseCov*diag(sqrt(RATES));
    else,
        % Return the noise covariance matrix
        % Useful in the moment equations, defined similarly to the
        % extended Kalman filter.
        noise = model.noiseCov*diag(RATES)*model.noiseCov';
    end

    % Safety flag
    %   0 = Do not test at all
    %   1 = Test and silently patch unexpected numerical situations
    %   2 = Test and warn on unexpected numerical situations
    %   3 = Test and Hard-fail on unexpected numerical situations
    %
    if model.safety>0,
        assertMatrixCondition(noise,model.safety,0);
    end


    '''#STOPCODE
