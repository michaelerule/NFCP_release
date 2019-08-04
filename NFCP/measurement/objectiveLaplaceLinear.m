function [f,g,Dh]=objectiveLaplaceLinear(...
    x,M,F,y,reg,bias,gain,variances)
    %{
    Objective function for a linear relationship between intensity
    :math:`\lambda` and latent
    state :math:`x`, with gain parameter :math:`m` and bias :math:`b`
    
    .. math::
        \lambda = (x+\text{bias})\cdot\text{gain}
    
    The negative log-likelihood (up to a constant)
    :math:`-\mathscr{L}` and its first two derivatives are:
    
    .. math::
        \begin{eqnarray}
        -\mathscr{L}  &=&
        \Delta_x (mx+b) - y \ln(mx+b) \\
        -\tfrac {\partial \mathscr{L} }{\partial x}  &=&
        \Delta_x m - y \left( \tfrac m {mx + b} \right)  \\
        -\tfrac {\partial^2 \mathscr{L} }{\partial x^2}  &=&
        y \left(\tfrac m {mx + b} \right)^{2}
        \end{eqnarray}
    
    In the SIRS model (QAR model in the case of neurons), the
    A (I) field is interpreted as the active state giving rise to
    point-process observations.
    
    This objetive jointly optimizes all variables (including those not
    directly observed), with a positivity constraint.
    
    Parameters
    ----------
    x : column vector length nspecies×n^2
        Guess for the vector of contentrations/intensities/rates
        at which to evaluate the objective function
    M : column vector length nspecies×n^2
        Prior mean for contentrations/intensities/rates
    F :
        Matrix such that prec = F'×F, can be low rank. Pre-factoring
        the precision matrix can accelerate objective computation.
    y  : column vector length nspecies×n^2
        Binned counts of point observations.
    reg :
        reg1/x regularization creating an asymptote at zero to
        forbid negative concentrations/intensities/rates.
        If this is nonzero, using a specialized gradient solver that
        avoids taking steps across zero is recommended.
    bias :
        A DC bias coupling x to point observations
    gain :
        A gain adjustment coupling x to point observations
        nonzero gain.
    variance : 
        Experimental: variance correction to expected likelihood.
        (set to zero);
    
    Returns
    -------
    f : scalar
        objective function (negative log probability)
    g : vector; optional
        gradient (only if needed)
    Dh : vector; optional
        hessian diagonal update; hessian = prec + diag(Dh)
        (only if needed)
    %}

    % ensure inputs are column vectors
    x = x(:);
    y = y(:);
    M = M(:);
    
    % get number of dimensions; m = n²
    m = numel(y);
    
    % apply gains/biases to latent states
    xb = x(m+1:2*m) + bias;
    xr = xb .* gain;
    
    % Only need the log-rate where `y` is nonzero (where spikes are 
    % observed). This can speed things up slightly. If using interpolated
    % (or otherwise non-square) basis functions, we might have fractional
    % values of spikes (<1) (a linear approximation). 
    % Truncate the smallest meaningful "spike" at 1e-1;
    sp = (y>1e-3)&(gain>0.0);
    
    % Expected log-likeliood coefficients
    c0      = 1./xb;
    c0(~sp) = 0.0;
    c1      = -c0.*c0;
    
    % Derivatives of log-likelihood
    lpyx0      = y.*log(xr);
    lpyx0(~sp) = 0.0;
    lpyx0      = lpyx0 - xr;
    lpyx2      = y.*c1;
    
    % second-order variance correction
    vc = 0.5.*variances;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Objective
    xM    = x-M ; % mean-centered x
    FxM   = F*xM; % use factorized precision matrix to compute prior
    x1    = 1./x; % regularization potential to prevent x<0
    f     = 0.5.*sum(FxM.^2)       ... % Prior
          + reg.*sum(x1)           ... % Regularization
          - sum(lpyx0 + vc.*lpyx2);    % Measurement

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gradient
    if nargout>1
        c2    = -2.*c1.*c0;
        lpyx1 = y.*c0 - gain;
        lpyx3 = y.*c2;
        pxM   = F'*FxM;
        x2    = x1.*x1;    % 2/x^2;
        rr    = x2.*-reg; % regularization potential to prevent x<0
        g     = pxM                 ... % Prior
              + rr                 ; % Regularization
        g(m+1:2*m) = g(m+1:2*m)      - (lpyx1 + vc.*lpyx3);    % Measurement
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Diagonal update to precision matrix needed to form Hessian
    % Hessian is H = prec + diag(update);
    if nargout>2
        c3    = -3.*c2.*c0;
        lpyx4 = y.*c3;
        x3    = x2.*x1;          % 1/x^3
        rr2   = 2.0.*(reg.*x3); % regularization potential to prevent x<0
        Dh    = rr2 ; % Regularization
        Dh(m+1:2*m) = Dh(m+1:2*m)  - (lpyx2 + vc.*lpyx4);    % Measurement
    end

