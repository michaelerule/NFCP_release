function [pxy,pM,pV,pP,dM,dP,dPM] = integrateMarginals(M,P,x,logpyx)
    %{
    In parallel, estimate moments of marginals by integrating
    
    Parameters
    ----------
    M : vector
        prior means
    P : vector
        prior precisions (1/variance)
    x : domain, e.g. linspace(0,1,100) if each marginal is in [0,1]. Must 
        be evenly spaces for probability mass to be counted properly. 
        Alternatively, x represents discrete (categorical) states or 
        probability bins.
    logpyx : Poisson spiking observation likelihoods
    
    Returns
    -------
    pxy : matrix
        log-posterior estimate evaluated at x, up to a constant, for all
        regions/channels/dimensions/variables
    pM : vector
        posterior mean estimates
    pV : vector
        posterior variance estimates
    pP : vector
        posterior precision estimates
    dM : vector
        change in mean
    dP : vector
        chang in precision
    dPM : vector
        premultiplied change in mean, precision
    %}
    
    minv = 1e-4; % minimum variance
    minp = 1e-7; % minimum probability per bin
    
    P = P(:); % ensure column vector
    M = M(:);
    x = x(:); % (:) ensure is nsamples x 1 column vector

    nvars    = length(M);
    nsamples = length(x);
    
    %V = max(V,minv);  % Ensure variances not pathologically small
    %P = 1./V;         % Precisions (inverse variances)
    P = max(P,minv);
    V = 1./P;

    % Add log-prior, log-likelihood to get log-posterior (up to constant)
    logpx = priorSurface(M,P,x);
    logpxy = logpyx + logpx;
    
    % Exponentiate and normalize to get posterior
    pxy = normrows(exp(logpxy));
        
    % Integrate to get posterior moments
    pM = sum(bsxfun(@times,pxy,x'),2);
    pV = sum(bsxfun(@times,pxy,(outer(x,pM,@minus)).^2'),2);
    
    % Un... truncate ?!? (approximately)
    Vadj = pV/(1-2/pi);
    Madj = pM - sqrt(Vadj)*sqrt(2/pi);
    Madj = max(Madj,0);
    %Vadj = sqrt(Vadj .* (pM  - Madj).^2*(pi/2));
    pV = Vadj;
    pM = Madj;
    
    % Also estimate the effect of the truncated domain on its own
    % this introduces a bias we may want to account for
    %pxx = normrows(exp(logpx));
    %bM = sum(bsxfun(@times,pxx,x'),2);
    %bV = sum(bsxfun(@times,pxx,(outer(x,bM,@minus)).^2'),2);
    %bV = max(bV,minv);
    %bP = 1./bV;
    
    assert(all(size(pV)==size(V)));
    
    % Ensure non-pathological variances and convert to precision
    pV = max(pV,minv);
    pP = 1./pV;
    
    % Ensure that all precisions increased.
    %assert(all(pP(:)-P(:)>-0.2));
    
    % For expectation-propagation, we need to know about the prior and
    % posterior precision -- but also need to keep track of a 
    % "divided-out" gaussian update.
    % Ppost = Pprior + Pdelta
    dP  = max(pP-P,1e-9);
    pP  = P + dP;
    dPM = pP.*pM - P.*M;% - bP.*bM;
    dM  = dPM./dP;
    
    
    
