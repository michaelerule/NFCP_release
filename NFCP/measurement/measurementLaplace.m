function [loglikelihood,pmode,pprec] = measurementLaplace(model,M,P,xypoints)
    %{
    Laplace approximation of measurement update.
    
    This routine explicitly constrains the result to be positive for ALL
    latent variables, not just those directly couplde to the point-process
    observations.    
    
    This update function is selected by setting `model.update` to 'Laplace'
    
    Parameters
    ----------
    model: struct
        Pre-initialized model struct. See `initializeModel`.
    M: vector
        Means for concentrations of species Q, A, R (concatenated)
    P: positive definite matrix
        Full precision (inverse covariance) matrix for the
        system. The first two states should be the Quiescent (Q) and
        actively firing (A) states. Each species has n x n spatial grid 
        packed in row-major order
    xypoints : npoints x 2 matrix
        Point xypoints in the form of a list of (x,y) locations.
        All points should like in the [0,1]² unit square.
        The function will also behave correctly if xypoints have been
        binned beforehand to an n x n grid.
        
    Returns
    -------
    loglikelihood : 
        The total log-likelihood of the xypoints after measurement 
        update, excluding constant factors, computed via Laplace
        approximation.
    pmode : 
        the posterior mode concentration for each location
    pprec : 
        the posterior precision for each location
    %}
    
    % Create a histogram of 2D xypoints (the observation counts)
    y = binCounts(model,xypoints);  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Newton-Raphson to find the posterior mode
    M    = M(:);  % Switch to column vector
    y    = y(:);  % Observations
    step = 1;     % Newton-Raphson step size
    f0   =-1/eps; % Stores previous best likelihood
    f    = 1/eps; % Stores current best likelihood
    % Clean up prior mean rates
    minr      = model.minrate;
    M(M<minr) = minr;      
    x         = M(:); % Estimated rates (initial guess)
    
    % Extract cholesky factor of prior precision matrix
    assertMatrixCondition(P,model.safety);
    F = chol(P); % x=chol(x)'*chol(x)

    % Get objective function
    ofun = getObjective(model.link);

    % Variance correction (not used presently)
    v = y.*0;

    gain = model.adjusted_gain;
    bias = model.adjusted_bias;
    objective = @(x) ofun(x,M,F,y,minr,bias,gain,v);

    % indices used to assign to diagonal covariance update (hessian); 
    [nRows,nCols] = size(P);
    pdiag = diag(P);
    H = P;

    % Apply Newton-Raphson method
    % Adaptive stepping:
    %   Values must be positive; negative values -> smaller step size
    %   Objective must monotonically increase, if not -> smaller step size
    %   Everything OK? -> try a larger step size
    iter = 1;

    while abs(f-f0)>abs(model.tol*f)   
        % Store old likelihood for comparison
        f0 = f;
        % Get Hessian, gradient, and addition to Hessian diagonal
        % for Newton-Raphson update
        [f,g,Dh] = objective(x); 
        H(1:(nRows+1):nRows*nCols) = pdiag + Dh;
        ch    = chol(H)   ; % x=ch'*ch
        delta = ch\(ch'\g);
        % Try to take a larger step size, then back-track if we
        % do not end up in a better place.
        step = step*20;   
        iterstep=1;      
        while true     
            %  Check that line search isn't stuck
            iterstep = iterstep+1;
            if (iterstep>model.maxiter*10), % too many iterations!
                if model.verbosity>2, fprintf(2,'Line search exceeded %d iterations\n',iterstep); end
                break;
            end
            if step<sqrt(eps), % line search failed
                if model.verbosity>2, fprintf(2,'Cannot improve objective along gradient\n'); end
                break;
            end
            % Check that the step size has not reduced to zero
            % assert_step(step,model.safety);
            % Move to next location
            newx = x-step*delta; 
            % Retry with smaller step if any intensities below minimum
            if any(newx<=model.minrate), step = step*0.5; continue; end
            % Check likelihood value and stop if we've improved
            f = objective(newx); 
            if f<=f0, break; end;
            % If not improving, try a smaller step
            step = step*0.5;
        end
        if any(newx<=model.minrate),
            model.error('Could not locate suitable non-negative concentration in update');
            newx = x;
        end
        if model.verbosity>2, display(['Current objective ' num2str(f)]); end;
        assert_vector_condition(newx,model.safety)
        x = newx; % Use this value for next iteration
        iter = iter+1;
        if (iter>model.maxiter), break; end
    end
    % Get Hessian, gradient, and addition to Hessian diagonal
    % for Newton-Raphson update. The hessian is used as the 
    % posterior precision matrix (inverse covariance)
    [f,g,Dh] = objective(x);
    H = P+diag(Dh);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Assign output arguments
    pmode = x;
    pprec = H;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimate of the likelihood
    loglikelihood = logLikelihood(model,M,P,pmode,pprec,y);
    
    assert(~isnan(loglikelihood));
    assert(all(isfinite(loglikelihood)));

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

function ofun=getObjective(link),
    %{ 
    Get objective function; 
    %}
    if strcmp(link,'linear'),
        ofun = @objectiveLaplaceLinear;
    elseif strcmp(link,'log'),
        % log-link is untested
        ofun = @objectiveLaplaceLog;
        fprintf(2,'Warning, the log-link function is untested');
    else
        error('Link function must be either linear or log');
    end











