function [logpyx] = likelihoodSurface(model,x,y)
    %{
    
    Parameters
    ----------
    model
    x
    y
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Log probability of observations, up to a constant, 
    % Model gains and biases enter into likelihood calculation

    gain = model.premultiplied_gain(:); 
    bias = model.adjusted_bias(:)./gain;
    bias(gain<1e-9) = 0.0;
    
    unobserved = (~isfinite(y)) | (gain<1e-9);
    y(unobserved) = 0.0;

    sp = (y>1e-3)&(gain>0.0)&(~unobserved);
    
    if strcmp(model.link,'linear'),
        % Linear Poisson log-likelihood with gain and bias
        % rate = gain*(bias+x)
        % (in log-probability the gain becomes an additive constant
        % and may be omitted)
        % log-rate contribution (applies when spikes observed)
        %rate           = bsxfun(@times,bsxfun(@plus,x',bias),gain);
        rate = bsxfun(@times,outer(bias,x,@plus),gain);
        assert(all(rate(:)>0));
        lxy            = bsxfun(@times,log(rate),y);
        lxy(~sp,1:end) = 0.0; % turn this off for no-spike regions
        % negative expected rate penalty (always applies)
        logpyx         = bsxfun(@minus,lxy,rate);
    elseif strcmp(model.link,'log'),
        % Log-Gaussian Poisson
        % rate = gain*(exp(x)+bias)
        % multiplicative correction of gain*exp(x)
        % in log transform becomes additive correction of +x
        rate = bsxfun(@times,outer(bias,exp(x),@plus),gain);
        assert(all(rate(:)>0));
        lxy            = bsxfun(@times,log(rate),y);
        lxy(~sp,1:end) = 0.0; % turn this off for no-spike regions
        % negative expected rate penalty (always applies)
        logpyx         = bsxfun(@minus,lxy,rate); 
        logpyx         = bsxfun(@plus,lxy,x); 
    elseif strcmp(model.link,'squared'),
        % Square-root-Gaussian Poisson
        % rate = gain*(x^2+bias)
        % multiplicative correction of gain*2*|x|
        % in log transform becomes additive correction of log(|x|)
        % which diverges at zero so be careful!
        rate = bsxfun(@times,outer(bias,x.^2,@plus),gain);
        assert(all(rate(:)>0));
        lxy            = bsxfun(@times,log(rate),y);
        lxy(~sp,1:end) = 0.0; % turn this off for no-spike regions
        % negative expected rate penalty (always applies)
        logpyx         = bsxfun(@minus,lxy,rate); 
        assert(all(x(:)>0));
        logpyx         = bsxfun(@plus,lxy,log(abs(x))); 
    else
        error('Supported link functions are linear, log, and (chi-)squared');
    end
    
    % Set unobserved channels to flat
    logpyx(unobserved,1:end) = 0.0;
    
    % Normalize to prevent underflow, 
    logpyx = min0rows(logpyx);
