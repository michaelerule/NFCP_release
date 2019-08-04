function logpx = priorSurface(M,P,x)
    %{
    Construct log-proabability densities for marginals.
    
    Matches moment of zero-truncated Gaussian to provided gaussian moments
    %}
    
    minv = 1e-3;
    minp = 1e-9/length(x); % minimum probability per bin

    M = M(:);
    P = P(:);
    V = 1./P;
    
    Vadj = V/(1-2/pi);
    Madj = M - sqrt(Vadj)*sqrt(2/pi);
    Madj = max(Madj,0);
    %Vadj = sqrt(Vadj .* (M  - Madj).^2*(pi/2));
    
    %Vadj = max(minv,Vadj);

    V = Vadj;
    M = Madj;
    P = 1./V;

    % Evaluate log-priors up to a constant
    % end result should be nvars x nsamples in shape
    xM2 = outer(M,x,@minus).^2;
    logpx = -0.5*(bsxfun(@minus, bsxfun(@times, xM2, P), log(P) ));
    
    % Normalize to prevent underflow, 
    logpx = min0rows(logpx);
    
    % Normalize in probability space, and pad out small probabilities
    % to avoid having a prior that is zero anywhere. the end result is
    % not normalized.
    % logpx is nvars x nsamples in shape
    px = normrows(exp(logpx));
    px = normrows(max(px,minp));
    logpx = log(px);
    
    % TODO: clean up
    % Normalize to prevent underflow, 
    logpx = min0rows(logpx);
end
