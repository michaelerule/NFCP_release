function [f,g,Dh]=objectiveSubspaceLinear(x,M,F,y,reg,bias,gain,variances)
    %{
    
    %}
    
    % ensure inputs are column vectors
    x = x(:);
    y = y(:);
    M = M(:);
    
    ex = exp(x);
    xb = ex  + bias;
    xr = xb .* gain;
    gx = gain .* ex;
    
    % Only need the log-rate where `y` is nonzero (where spikes are 
    % observed). This can speed things up slightly. If using interpolated
    % (or otherwise non-square) basis functions, we might have fractional
    % values of spikes (<1) (a linear approximation). 
    % Truncate the smallest meaningful "spike" at 1e-3;
    sp    = (y>1e-3)&(gain>0.0);
    
    % Expected log-likeliood coefficients
    c0 = ex./xb;
    c0(~sp) = 0.0;
    c1 = c0.*(1-c0);
    
    % Derivatives of log-likelihood
    lpyx0 = y.*log(xr);
    lpyx0(~sp) = 0.0;
    lpyx0 = lpyx0 - xr;
    lpyx2 = y.*c1 - gx;
    
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
        c2    = c1.*(1-2*c0);
        lpyx1 = y.*c0 - gx;
        lpyx3 = y.*c2 - gx;
        pxM   = F'*FxM;
        x2    = x1.*x1;    % 2/x^2;
        rr    = x2.*-reg; % regularization potential to prevent x<0
        g     = pxM                 ... % Prior
              + rr                  ... % Regularization
              - (lpyx1 + vc.*lpyx3);    % Measurement
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Diagonal update to precision matrix needed to form Hessian
    % Hessian is H = prec + diag(update);
    if nargout>2
        c3    = c1.*(1-6*c1);
        lpyx4 = y.*c3 - gx;
        x3    = x2.*x1;          % 1/x^3
        rr2   = 2.0.*(reg.*x3); % regularization potential to prevent x<0
        Dh    = rr2                 ... % Regularization
              - (lpyx2 + vc.*lpyx4);    % Measurement
    end
    
    
    
