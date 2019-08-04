function [pMa,pPa] = laplaceUpdateNewtonRaphson(model,Ma,Pa,Ca,chPa,y)
    %{
    Solve the measurement update using the Laplace approximation using
    a customized Newton-Raphson algorithm
    
    Note: verified that bias and gain parameterization matches that in
    the logLikelihood function 25 July 2018
    
    Parameters
    ----------
    model: struct
        Pre-initialized model struct. See `initializeModel`.
    Ma: vector
        Prior means for active state
    Pa: positive definite matrix
        Prior precision (inverse covariance) matrix for active state 
    Ca: positive definite matrix
        Prior covariance matrix for active state 
    y : vector
        Binned and adjusted count observations
        
    Returns
    -------
    pMa : vector
        Posterior means
    pPa : positive definite matrix
        Posterior precisions
    %}
    
    m    = model.nn; % No. basis function sper state
    step = 1;        % Newton-Raphson step size
    f0   =-1/eps;    % Stores previous best likelihood
    f    = 1/eps;    % Stores current  best likelihood
        
    % get objective function
    ofun = getObjective(model);
    
    v = diag(Ca).*0;
    use_expected_log_likelihood = false;

    % Clip prior rates to minimum rate
    minr = model.minrate; 
    Ma = max(Ma,minr);
    
    gain = model.adjusted_gain;
    bias = model.adjusted_bias;
    objective = @(x) ofun(x,Ma,chPa,y,minr,bias,gain,v);
    
    % Initial guess at prior mean
    x = Ma;
    
    % Estimate likelihood curvature without positivity barrier (minr=0),
    % and with zero variance. This is the "correct" curvature at the
    % estimated posterior mode.
    v0 = v.*0;
    curvature = @(x) ofun(x,Ma,chPa,y,minr,bias,gain,v0);

    % indices used to assign to diagonal covariance update (hessian); 
    [nRows,nCols] = size(Pa);
    pdiag         = diag(Pa);
    hessian       = Pa + eye(m).*model.reg_inverse;
    H2            = hessian;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Newton-Raphson
    % Adaptive stepping:
    %   Values must be positive; negative values -> smaller step size
    %   Objective must monotonically increase, if not -> smaller step size
    %   Everything OK? -> try a larger step size
    iter = 1;
    while abs(f-f0)>abs(model.tol*f)   
        % Store old likelihood for comparison
        f0 = f;
        
        % Get Hessian, gradient, and addition to Hessian diagonal
        [f,g,Dh] = objective(x); 
        %Dh = min(max(Dh,0),1e12);
        hessian(1:(nRows+1):nRows*nCols) = pdiag + max(Dh,0);
        %{
        if rcond(hessian)<1e-15,
            % This seems to be caused by very large curvatures
            % (not very small ones). 
            error('hessian is ill-conditioned');
            while rcond(hessian)<1e-15,
                hessian = hessian + eye(m).*1e-3;
            end
        end
        %}
        ch    = chol(hessian); % x=ch'*ch
        delta = ch\(ch'\g);
        
        %{
        if use_expected_log_likelihood,
            % Use Laplace estimated covariance to update the posterior 
            % variance, which is used to approximate the expected ll
            % this is not guaranteed to converge
            % This should be combined with a Dkl penalty to approximate
            % a variational update.
            [f2,g2,Dh2] = curvature(x);
            H2(1:(nRows+1):nRows*nCols) = pdiag + max(Dh2,0);
            chCest = (chol(H2)\model.Inn)';
            v      = diag(chCest'*chCest);
            objective  = @(x) ofun(x,Ma,chPa,y,minr,ab,ag,v);
        end
        %}
               
        % Try to take a large step size, back-track if it gets worse
        step     = step*20;   
        iterstep = 1;      
        while true     
            %  Check that line search isn't stuck, and abort if it is.
            iterstep = iterstep+1;
            if (iterstep>model.maxiter*10), % Something is wrong
                if model.verbosity>2, 
                    fprintf(2,'Line search exceeded %d iterations\n',iterstep);  
                end
                break;
            end
            if step<sqrt(eps), % Something is wrong
                if model.verbosity>2, 
                    fprintf(2,'Cannot improve objective along gradient\n');  
                end
                break;
            end
            % Check that the step size has not reduced to 0 (bad gradient)
            assert_step(step,model.safety);
            newx = x-step*delta;   % Move to next location
            % Retry with smaller step if any intensities fall below minimum
            if any(newx<=minr), step = step*0.5; continue; end
            f = objective(newx);   % Check likelihood, and
            if f<=f0, break; end;  % stop if we've improved,
            step = step*0.5;       % If not improving, try a smaller step
        end 
        % In rare cases line-search can fail while searching for a 
        % non-negative concentration in the updates
        if any(newx<=minr),
            model.error('Could not find positive intensity update');
            newx = x;
        end
        if model.verbosity>2, display(['ll~' num2str(f)]);  end;
        assert_vector_condition(newx,model.safety)
        x    = newx;   % Use this value...
        iter = iter+1; % for the next iteration
        if (iter>model.maxiter), break; end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Unpacking
    % Get Hessian, gradient, and addition to Hessian diagonal
    % The hessian is used as the posterior precision matrix
    [f,g,Dh] = curvature(x);
    %[f,g,Dh] = objective(x);
    hessian(1:(nRows+1):nRows*nCols) = pdiag + max(Dh,0);
    pMa   = x;
    pPa   = hessian;
    
    pCa   = cinv(pPa);
    pVa   = diag(pCa);
    pSa   = pVa.^0.5;
    
    % Use the curvature at the estimated posterior mode as an 
    % approximation of the inverse covariance (precision)
    % For the Poisson log-likelihood, the curvature is very steep
    % near zero, and the distribution can be very skewed. In this case
    % the local curvature can severely under-estimate the variance
    % we make a small adjustment to the location of the posterior mode
    % in attempt to address this issue, considering that on average 
    % most of the probability mass will be to the right of the mode.
    % This adjustment uses pSa, which is the posterior standard deviation.
    % DISABELING FOR NOW, VARIANCE IS TOO HIGH 
    x_adjusted = x;% + pSa.*sqrt(2/pi);
    [f,g,Dh] = curvature(x_adjusted);
    %[f,g,Dh] = objective(x);
    hessian(1:(nRows+1):nRows*nCols) = pdiag + max(Dh,0);
    pMa   = x;
    pPa   = hessian;
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error handlers

function assert_step(step,safety,tolerance)
    %{
    If the delta vector does not point down-gradient, no step size
    will improve the objective function; This is probably because 
    the gradient and/or hessian are wrong. 
    
    Parameters
    ----------
    step : scalar
        Minimum newton-raphson step size above which no improvement
        in objective can be found; if this gets to zero within
        machine precision, it means that the gradient is not valid. 
    safety : int; default 2
        Safety flag
          2 = Test and warn on unexpected numerical situations
          3 = Test and Hard-fail on unexpected numerical situations
    tolerance : scalar; optional, default sqrt(eps)
        minimum ratio between smallest and largest magnitude to be
        allowed
    %}
    if nargin<3, tolerance=sqrt(eps); end;
    if (step<=eps)
        msg = sprintf([...
            'Measurement update: unable to improve objective\n'...
            'function for Laplace approximation in the\n'...
            'direction of the gradient; something is wrong']);
        if safety==3,
            error(msg);
        elseif safety==2,
            warning(msg);
        end
    end

function assert_vector_condition(x,safety,tolerance)
    %{
    A large range in x will cause the Hessian to become singular
    up to machine precision; this must be avoided. Raises an error
    if the ratio of smallest to largest value is less than the 
    tolerance.
    
    Parameters
    ----------
    x : vector
        absolute magnitudes will be checked
    safety : int; default 2
        Safety flag
          2 = Test and warn on unexpected numerical situations
          3 = Test and Hard-fail on unexpected numerical situations
    tolerance : scalar; optional, default eps
        minimum ratio between smallest and largest magnitude to be
        allowed
    %}
    if nargin<3, tolerance=eps; end;
    if min(abs(x))/max(abs(x))<tolerance
        msg = sprintf([...
            'Measurement update: exceeded %0.2e size difference\n'...
            'between smallest and largest intensity; This will\n'...
            'lead to a singular covariance model. Aborting.\n'...
            'Minimum is %e, maximum is %e. To control large\n'...
            'values, increase L2 regularization'],...
            tolerance,min(abs(x)),max(abs(x))  );
        if safety==3,
            error(msg);
        elseif safety==2,
            warning(msg);
        end
    end  


function ofun = getObjective(model),
    %{
    get objective function; returns likelihood, gradient, hessian
    %}
    if strcmp(model.link,'linear'),
        %ofun = @objectiveSubspaceLinear;
        ofun = @objectiveSubspaceLinear;
    elseif strcmp(model.link,'log'),
        ofun = @objectiveSubspaceLog;
        fprintf(2,'The log link function is not tested\n');
    elseif strcmp(model.link,'squared'),
        ofun = @objectiveSubspaceSquared;
        fprintf(2,'The squared link function is not tested\n');
    else
        error('Supported link functions for the Laplace update on the active subspace are linear, log, and (chi-) squared');
    end


