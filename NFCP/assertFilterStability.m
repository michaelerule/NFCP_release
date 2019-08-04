function assertFilterStability(model,M,P,error_mode)
    %{
    Determine whether the state and measurement update predict rate
    increases even when there are no observed spikes.
    
    Parameters
    ----------
    model: ``struct``
        Pre-initialized model struct. See `initializeModel`.
    M: ``vector``
        Means for intensities of states. The first states must represent
        Quiescent (Q) neurons. The second states must represent actively
        firing (A) neurons. The remaining states can implement arbitrary
        linear dynamics (with Poisson noise). 
    P: ``matrix``
        Precision matrix
    error_mode: ``int`` ∈[0,1,2]
        Error mode flag default is 0
        0 : print warning to standard error, do not throw warning
        1 : raise warning (will appear with stack trace)
        2 : raise error
        def
    %}
    
    % Default to printing warning messages to stderr
    if nargin<4,
        error_mode = 0;
    end
    
    % Generate vector representing observation model (no bias)
    m = model.nn;
    O = zeros(model.nn,1);
    M = [O; O+model.gain*model.dx; O];
    
    % Compute closed-form update in absence of spikes
    ch = chol(P);
    delta = - ch\(ch'\M);
    M_closedform = M(:) + delta(:);
    
    % Check mean-update to get the rate of firing increase
    M_meanupdate = meanUpdate(model,M_closedform,P) + delta(:);

    % Verify that estimated rates do not increase in absense of
    % spiking observations.
    A_meanupdate   = M_meanupdate(m+1:2*m);
    A              = M(m+1:2*m);
    maxincrease    = max(A_meanupdate(:)-A(:));
    
    % Reducing excitation, increasing the rate at which bursting cells
    % become refractory, and increasing either the intrinsic or 
    % regularizing noise, can all help to stabilise a filter that is
    % otherwise unstable (i.e. predicting firing rate increases when
    % spikes are not observed)
    if maxincrease>0,
        message = 'Large self-excitation, reduce rAA to stabilize filter?';
        if error_mode==0,
            fprintf(2,'%s\n',message);
        elseif error_mode==1,
            warning(message);
        elseif error_mode==2,
            error(message);
        end
    end
    
    
    
