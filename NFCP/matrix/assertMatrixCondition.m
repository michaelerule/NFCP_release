function assert_matrix_condition(H,safety,tolerance)
    %{
    If a precision matrix is near-singular, we cannot find a gradient.
    This checks the reciprocal condition number of the precision
    or Hessian matrix. If it is smaller than the tolerance, an error
    is raised.
    
    Parameters
    ----------
    H : matrix
        Precision or Hessian matrix from Newton-Raphson update
        for Poisson likelihood
    safety : int; default 2
        Safety flag
          0 = Do nothing (skip all checks)
          1 = Repair errors when possible (this function does not repair)
          2 = Test and warn on unexpected numerical situations
          3 = Test and Hard-fail on unexpected numerical situations
    tolerance : scalar; optional, default `eps`
        Minimum reciprocal condition number permitted. If rcond(H)
        is smaller than this, an error will be thrown. Defaults to
        `eps`.
    %}
    if nargin<2, 
        safety=2; 
    end;
    if nargin<3, 
        tolerance=eps; 
    end;
    
    % Skip checks for low safety levels
    if safety>=2,
        if ~isreal(H),
            msg = 'Unexpected complex-valued matrix encountered';
            if safety==2,
                warning(msg);
            elseif safety==3,
                error(msg);
            end
        end
        if rcond(H)<tolerance
            msg = sprintf('Non-positive matrix rcond %e < %e',rcond(H),tolerance);
            if safety==2,
                warning(msg);
            elseif safety==3,
                error(msg);
            end
        end
    end
