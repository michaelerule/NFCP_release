function [f,g,Dh]=objectiveSubspaceSquared(x,M,F,y,reg,bias,gain,variances)
    %{
    
    
    
    %}
    
    % enable assertions
    debugmode = true;
    
    % ensure inputs are column vectors
    x = x(:);
    y = y(:);
    M = M(:);
    
    % interpret negative the same as positive, since value is to be
    % squared
    %x = abs(real(x));
    
    ex = x.*x;
    xb = ex  + bias;
    xr = xb .* gain;
    gx = gain .* ex;
    
    if debugmode,
        assert(all(xb>=0));
        assert(all(xr>=0));
        assert(all(gx>=0));
    end
    
    % Only need the log-rate where `y` is nonzero (where spikes are 
    % observed). This can speed things up slightly. If using interpolated
    % (or otherwise non-square) basis functions, we might have fractional
    % values of spikes (<1) (a linear approximation). 
    % Truncate the smallest meaningful "spike" at 1e-3;
    sp    = (y>1e-3)&(gain>0.0);
    
    % Expected log-likeliood coefficients
    d               = 2./xb;
    d(~sp)          = 0.0;
    d(~isfinite(d)) = 0.0;
    c0              = d.*x;
    c1              = d-c0.^2;
    
    % Derivatives of log-likelihood
    lpyx0 = y.*log(xr);
    lpyx0(~sp) = 0.0;
    lpyx0(~isfinite(lpyx0)) = 0.0;
    lpyx0 = lpyx0 - xr;
    lpyx2 = y.*c1-2*gain;
    
    if debugmode,
        assert(all(isfinite(d)));
        assert(all(isfinite(c0)));
        assert(all(isfinite(c1)));
        %assert(all(c1<=0.0)); % convexity condition
        assert(all(isfinite(lpyx0)));
        assert(all(isfinite(lpyx2)));
    end
    
    % second-order variance correction
    vc = 0.5.*variances;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Objective
    xM    = x-M ; % mean-centered x
    FxM   = F*xM; % use factorized precision matrix to compute prior
    f     = 0.5.*sum(FxM.^2)       ... % Prior
          - sum(lpyx0 + vc.*lpyx2);    % Measurement
    if reg>0,
        x1 = 1./x;              % regularization potential to prevent x<0
        x1(~isfinite(x1)) = 0.0;
        f  = f + reg.*sum(x1);  % Regularization
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gradient
    if nargout>1
        c2    = d.*c0.*(2.*x.*c0-3);
        lpyx1 = y.*c0-2*gain.*x;
        lpyx3 = y.*c2;
        pxM   = F'*FxM;
        g     = pxM              ... % Prior
              - (lpyx1 + vc.*lpyx3); % Measurement
        if reg>0,
            x2 = x1.*x1;   % 2/x^2;
            rr = x2.*-reg; % regularization potential to prevent x<0
            g  = g + rr;   % Regularization
        end
        if debugmode,
            assert(all(isfinite(c2)));
            assert(all(isfinite(lpyx1)));
            assert(all(isfinite(lpyx3)));
            assert(all(isfinite(g)));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Diagonal update to precision matrix needed to form Hessian
    % Hessian is H = prec + diag(update);
    if nargout>2
        c3    = 3*d.^2 - 6*c1.^2;
        lpyx4 = y.*c3;
        Dh    = -(lpyx2+vc.*lpyx4);% Measurement
        if reg>0,
            x3  = x2.*x1;          % 1/x^3
            rr2 = 2.0.*(reg.*x3);  % regularization potential to prevent x<0
            Dh  = Dh + rr2;        % Regularization
        end
        if debugmode,
            assert(all(isfinite(c3)));
            assert(all(isfinite(lpyx4)));
            assert(all(isfinite(Dh)));
            %assert(all(Dh>=0.0));       % convextiy
        end
end
