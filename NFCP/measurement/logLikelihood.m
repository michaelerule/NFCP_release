function [loglikelihood] = logLikelihood(model,M,P,x,pP,y)
    %{
    Note: verified that bias and gain parameterization matches that in
    the laplaceUpdateNewtonRaphson function 25 July 2018
    
    Note: y cound values should be pre-multiplied by model.alpha
    
    Parameters
    ----------
    model : model struct
    M : prior means
    P : prior precisions
    x : posterior means
    pP: posterior precisions
    y : count observations
    
    Returns
    -------
    loglikelihood : scalar or vector depending on settings
    %}
    
    
    %{
    extract gains and biases
    
    +--------------------------+-------------------------------------------+
    | Variable                 | Units                                     | 
    +==========================+===========================================+
    | model.adjusted_bias      | unitless                                  |
    +--------------------------+-------------------------------------------+
    | model.adjusted_gain      | alpha∙spikes /   bin_x²∙bin_t             |
    +--------------------------+-------------------------------------------+
    %}
    gain = model.adjusted_gain(:);
    bias = model.adjusted_bias(:);
    bias(gain<1e-9) = 0.0;
    
    % detect regions that are not observed
    unobserved = (~isfinite(y)) | (gain<1e-9);
    y(unobserved) = 0.0;
    sp = (y>1e-3)&(gain>0.0)&(~unobserved);
    
    % Full means/precisions if not using subspace update
    %  M : prior means
    %  P : prior precisions
    %  x : posterior means
    %  pP: posterior precisions
    priorMeanFull = M(:);
    postMeanFull  = x(:);
    priorPrecFull = P ;
    postPrecFull  = pP;
    
    % Extract active subspace that couples to point-process observations
    %  M : prior means
    %  P : prior precisions
    %  x : posterior means
    %  pP: posterior precisions
    if size(M(:),1)==model.nn,
        % Only A subspace passed to likelihood
        priorMeanA = M(:);
        postMeanA  = x(:);
        priorPrecA = P ;
        postPrecA  = pP;
    else,
        % Full model passed to likelihood, extract A subspace
        priorMeanA = model.getA'*M(:);
        postMeanA  = model.getA'*x(:);
        priorPrecA = model.getA'*P *model.getA;
        postPrecA  = model.getA'*pP*model.getA;
    end    
    
    chPpostA   = chol(postPrecA) ; % cholesky of posterior precision matrix
    chPpriorA  = chol(priorPrecA); % cholesky of prior precision matrix
    
    % Compute rates and likelihood curvature under various possible
    % link-functions. 
    % This computes a curvature correction based on the likelihood
    % curvature at the posterior mean.
    %
    % xb: the posterior rates, bias adjusted
    % xr: the posterior rates, bias and gain adjusted
    % lpyx2Post : the curvature at the posterior mode
    % xrPrior : the bias/gain adjusted rates at the PROIR mean
    if strcmp(model.link,'linear'),
        % Calculate for the prior mean
        xbPrior           = priorMeanA + bias;
        xrPrior           = xbPrior .* gain;
        c0Prior           = 1./xbPrior;
        c0Prior(~sp|~isfinite(c0Prior)) = 0.0;
        c1Prior           = -c0Prior.*c0Prior;
        lpyx2Prior        = y.*c1Prior;
        % Calculate for the posterior mean
        xbPost            = postMeanA + bias;
        xrPost            = xbPost .* gain;
        c0Post            = 1./xbPost;
        c0Post(~sp|~isfinite(c0Post)) = 0.0;
        c1Post            = -c0Post.*c0Post;
        lpyx2Post         = y.*c1Post;
    elseif strcmp(model.link,'log'),
        % Calculate for the prior mean
        exPrior           = exp(priorMeanA);
        xbPrior           = exPrior  + bias;
        xrPrior           = xbPrior .* gain;
        gxPrior           = gain .* exPrior;
        c0Prior           = exPrior./xbPrior;
        c0Prior(~sp|~isfinite(c0Prior)) = 0.0;
        c1Prior           = c0Prior.*(1-c0Prior);
        lpyx2Prior        = y.*c1Prior - gxPrior;
        % Calculate for the posterior mean
        exPost            = exp(postMeanA);
        xbPost            = exPost + bias;
        xrPost            = xbPost .* gain;
        gxPost            = gain .* exPost;
        c0Post            = exPost./xbPost;
        c0Post(~sp|~isfinite(c0Post)) = 0.0;
        c1Post            = c0Post.*(1-c0Post);
        lpyx2Post         = y.*c1Post - gxPost;
    elseif strcmp(model.link,'squared'),
        % Calculate for the prior mean
        exPrior           = priorMeanA.*priorMeanA;
        xbPrior           = exPrior  + bias;
        xrPrior           = xbPrior .* gain;
        gxPrior           = gain .* exPrior;
        dPrior            = 2./xbPrior;
        dPrior(~sp|~isfinite(dPrior))   = 0.0;
        c0Prior           = dPrior.*x;
        c1Prior           = dPrior-c0Prior.^2;
        lpyx2Prior        = y.*c1Prior-2*gain;
        % Calculate for the posterior mean
        exPost            = postMeanA.*postMeanA;
        xbPost            = exPost + bias;
        xrPost            = xbPost .* gain;
        gxPost            = gain .* exPost;
        dPost             = 2./xbPost;
        dPost(~sp|~isfinite(dPost)) = 0.0;
        c0Post            = dPost.*x;
        c1Post            = dPost-c0Post.^2;
        lpyx2Post         = y.*c1Post-2*gain;
    else
        error('Supported link functions are linear, log, and squared');
    end
    
    %{
    Possible likelihood approximations combine various approximations 
    of the data likelihood, and any penalties for inaccurate state 
    predictions reflected in a large change between the prior state
    estimate and the posterior state estimate after measurement. 
    %}
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Expected log-likelihood near posterior mean via 2nd-order expansion
    % Compute point likelihood estimate at the POSTERIOR mean. 
    % This estimate does not incorporate a covariance correction.
    lpyx0Post = y.*log(max(xrPost,1e-9));
    lpyx0Post(~sp|~isfinite(lpyx0Post)) = 0;
    lpyx0Post = lpyx0Post - xrPost;
    lpyx2Post(~isfinite(lpyx2Post)) = 0;
    
    % Compute covariance correction at the POSTERIOR mean
    % First, extract the diagonal (marginal) variances by inverting the
    % posterior precision matirx:
    chpCPost = (chPpostA\eye(size(chPpostA,1)))';
    pCPost   = chpCPost'*chpCPost;
    vcPost   = 0.5.*diag(pCPost);
    % Then, multiply this variance correcion by the likelihood curvature
    % at the posterior mean (and remove unphysical values)
    correctionPost = vcPost.*lpyx2Post;
    correctionPost(~isfinite(correctionPost)) = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Expected log-likelihood near prior mean via 2nd-order expansion
    % Compute point likelihood estimate at the prior mean. 
    % This estimate does not incorporate a covariance correction.
    lpyx0Prior = y.*log(max(xrPrior,1e-9));
    lpyx0Prior(~sp|~isfinite(lpyx0Prior)) = 0;
    lpyx0Prior = lpyx0Prior - xrPrior;
    lpyx2Prior(~isfinite(lpyx2Prior)) = 0;
    
    % Compute covariance correction at the PRIOR mean
    % First, extract the diagonal (marginal) variances by inverting the
    % posterior precision matirx:
    chpCPrior = (chPpriorA\eye(size(chPpriorA,1)))';
    pCPrior   = chpCPrior'*chpCPrior;
    vcPrior   = 0.5.*diag(pCPrior);
    % Then, multiply this variance correcion by the likelihood curvature
    % at the posterior mean (and remove unphysical values)
    correctionPrior = vcPrior.*lpyx2Prior;
    correctionPrior(~isfinite(correctionPrior)) = 0;
    %% ^^^ correctionPrior can sometimes fail ^^^
    % This occurs when the prior variance is too large, since the 2nd-
    % order approximation is only valid for small variance. 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    Also compute quantities based on the pre-update and
    post-update latent state distributions, which includes states
    that might not be directly observed (and therefore omitted from 
    the terms that depend on the data likelihood)
    %}
    D              = numel(x);
    chPpostFull    = chol(postPrecFull) ; % cholesky of posterior precision matrix
    chPpriorFull   = chol(priorPrecFull); % cholesky of prior precision matrix
    PriorPostFull  = chPpostFull \(chPpostFull' \postPrecFull ); % pP\P=C\pC;
    PostPriorFull  = chPpriorFull\(chPpriorFull'\priorPrecFull); % P\pP=pC\C;
    tracePriorPost = trace(PriorPostFull);
    tracePostPrior = trace(PostPriorFull);
    dmu            = postMeanFull - priorMeanFull;
    priorQuadratic = dmu'*priorPrecFull*dmu;
    postQuadratic  = dmu'*postPrecFull*dmu;

    % Log-determinants
    % The negative sum of the log-diagonal of the cholesky factor of a 
    % precision (inverse covariance) matrix is equal to one-half the
    % log-determinant of the corresponding covariance matrix
    % log-determinant of posterior COVARIANCE matrix
    logdetPost  = -2.*sum(log(max(model.llreg,diag(chPpostFull))));
    % log-determinant of prior COVARIANCE matrix
    logdetPrior = -2.*sum(log(max(model.llreg,diag(chPpriorFull))));
    logdetRatio = logdetPrior - logdetPost;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Likelihood via KL divergence of P from Q (ELBO)
    % log Py >= <log Py|a>_Q - Dkl(Q||Pr(A))
    % D_{KL}(Q || Pr(A)) can be computed exactly 
    DklQPa  = 0.5.*(logdetRatio - D + tracePriorPost + priorQuadratic);
    
    % <log Pr(Y|A) >_Q(a) can be computed via second-order approximation
    % using the previously calculated likelihood and curvature at the
    % POSTERIOR mean.
    logPyaQ = lpyx0Post + correctionPost;
    logPyaQ(~isfinite(logPyaQ)) = 0.0;
    logPyaQ = sum(logPyaQ);
    % Evidence lower bound is the expected log-likelihood under posterior Q
    % minus KL divergence of prior Pr(A) from Q
    ELBO = logPyaQ - DklQPa;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Expected log-likelihood approach (ELL)
    % log Py ~ <logPy|a>_Pa + Dkl(Pr(A)||Q)
    % D_{KL}(Pr(A) || Q) can be computed exactly 
    DklPaQ  = 0.5.*(-logdetRatio - D + tracePostPrior + postQuadratic);
    % <log Pr(Y|A) >_P(a) can be computed via second-order approximation
    % using the previously calculated likelihood and curvature at the
    % PRIOR mean.
    logPyaPa = lpyx0Prior + correctionPrior;
    logPyaPa(~isfinite(logPyaPa)) = 0.0;
    logPyaPa = sum(logPyaPa);
    % Expected log-likelihood under Prior Pr(A) plus KL divergence of prior
    % Q from Pr(A)
    ELL = logPyaPa + DklPaQ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Laplace approximated likelihoods (at prior and posterior mean)
    % Compute posterior Laplace-approximated liklihood
    lgPyLaplaceAtPostprior = -0.5.*(logdetRatio + priorQuadratic);
    lgPyLaplaceAtPost = lgPyLaplaceAtPostprior + sum(lpyx0Post);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Return all variants of likelihood estimate
    loglikelihood = [ELBO,ELL,lgPyLaplaceAtPost]';
    
    if strcmp(model.likemethod,'ML'),
        loglikelihood = [logPyaQ]';
    elseif strcmp(model.likemethod,'Laplace'),
        loglikelihood = [lgPyLaplaceAtPost]';
    elseif strcmp(model.likemethod,'ELBO'),
        loglikelihood = [ELBO]';
    elseif strcmp(model.likemethod,'ELL'),
        loglikelihood = [ELL]';
    elseif strcmp(model.likemethod,'debug'),
        % Components of variational likelihood separated
        loglikelihood = [logPyaQ,-DklQPa]';
    else
        error('likelihood method should be laplace, ELBO, ELL, or debug');
    end
       
    assert(all(~isnan(loglikelihood)));
    assert(all(isfinite(loglikelihood)));
    assert(all(isreal(loglikelihood)));
    
    
