function [f,g,Dh]=objectiveLaplaceLog(...
    x,M,F,y,reg1,reg2,bias,gain,volume)
    %{
    Considering a linear relationship between intensity
    :math:`\lambda` and latent
    state :math:`x`, with gain parameter :math:`m` and bias :math:`b`
    
    .. math::
        \lambda = xm+b
    
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
    A (I) field takes the role of latent variable x.
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
    reg1 :
        reg1/x regularization creating an asymptote at zero to
        forbid negative concentrations/intensities/rates.
        If this is nonzero, using a specialized gradient solver that
        avoids taking steps across zero is recommended.
    reg2 :
        L2 regularization on concentrations/intensities/rates
    bias :
        A DC bias coupling x to point observations
        TODO can be spatially inhomogeneous
    gain :
        A gain adjustment coupling x to point observations
        TODO can be spatially inhomogeneous
        TODO make measurement update only operate on variables with
        nonzero gain.
    volume : 
        alpha*dt*dx*dx
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Argument checking and precoMtation

    % ensure inputs are column vectors
    x = x(:);
    y = y(:);
    M = M(:);
    % get number of dimensions; m = n²
    m = numel(y);
    % rates for variables coupled to point observations
    
    lambda = bias.*exp(gain.*x(m+1:2*m));
    rate   = volume.*lambda;
    
    % We only need the log-rate where `y` is nonzero (where spikes are 
    % observed). This can speed things up slightly. If using interpolated
    % (or otherwise non-square) basis functions, we might have fractional
    % values of spikes (<1) (a linear approximation). 
    % Truncate the smallest meaningful "spike" at 1e-1;
    goodchannels = gain>0.0;
    sp = (y>0.1)&(goodchannels);
    logsp = sum(y(sp).*log(rate(sp)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Objective: most frequently computed value
    xM  = x-M ; % mean-centered x
    FxM = F*xM; % use factorized precision matrix to compute prior
    x1  = 1./x;
    f = 0.5 .*sum(FxM.^2) ...            % Prior term
      + sum(rate) - logsp ...        % Measurement update term
      + reg1.*sum(x1) + reg2.*sum(x.^2); % Regularization terms
    f = double(f);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gradient: compute only if necessary
    if nargout>1
        dd(~goodchannels) = 0.0;
        pxM = F'*FxM;
        x2 = x1./x; % x.^-2;
        g = pxM ...                     % Prior
          + (-reg1).*x2 + (reg2*2).*x;  % Regularization
        g(m+1:2*m) = g(m+1:2*m) + gain.*(rate-y); % Measurement
        g = double(g);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Diagonal update to precision matrix needed to
    % form the Hessian: compute only if necessary
    % To recover hessian compute
    % H = prec + diag(update);
    % (it may not be necessary to explicitly construct the above
    % matrix; or faster constructions involving in-place modification
    % of a copy of prec may be preffered; this is a tight-loop, so
    % memory optimization matter).
    if nargout>2
        x3 = x2./x;%x.^-3
        Dh = 2.0.*(reg1.*x3 + reg2); % Regularization
        Dh(m+1:2*m) = Dh(m+1:2*m)+(gain.*gain.*rate); % Measurement
        Dh = double(Dh);
end
