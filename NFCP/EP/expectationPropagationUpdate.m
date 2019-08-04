function [pMa,pPa] = expectationPropagationUpdate(model, Ma, Pa, Ca, y)
    %{
    Pseudocode 
    1. Get marginal priors
    2. Get marginal likelihoods (should be independent)
    3. Perform measurement on marginals by integration
    4. Compute divided-out change 
    5. Apply all divided-out changes jointly to multivariate priors
    
    Parameters
    ----------
    model : the model structure
    Ma : prior means
    Pa : prior precision
    Ca : prior covariances
    y  : spike counts, or NaN if no observations
    %}
    
    minv    = 1e-3;
    NEPITER = 150;
    if isfield(model,'maxiter'), NEPITER=model.maxiter; end;
    
    tol = 1e-9;
    Ta  = diag(Pa);           % Prior diagonal precision, \tau_a
    Va  = max(diag(Ca),minv); % Prior diagonal variance, \sigma^2_a

    % latent variable integration domain
    % needs to be specified somehow
    x = linspace(1e-9,max(Ma+10*sqrt(Va)),200);
    
    m = length(Ma);
    I = eye(m);
    
    logpyx = likelihoodSurface(model,x,y);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Expectation-Propagation (EP) iterations
    V = max(Va(:),minv);
    M = Ma(:);
    T = 1./V;
    for iteration=1:NEPITER, 
        % Compute Guassian-like updates (dM, dP, dMP) based on updates
        % computed by integrating marginals independently. 
        [pxy,pM,pV,pP,dM,dP,dPM] = integrateMarginals(M,T,x,logpyx);
        
        % Apply updates to joint prior distributions to estimate a joint
        % posterior
        pPa  = Pa + diag(dP);
        pPMa = Pa*Ma + dPM;
        pMa  = pPa\pPMa;
        
        % Convert back to covariance in order to get marginal variances
        chpP = chol(pPa); % pP = chpP'*chpP;
        chpC = (chpP\I)';
        pC   = chpC'*chpC;
        % Get marginal variances;
        pV   = diag(pC);
        pV   = max(pV,minv);

        % pT now contains propagated information from _all_ measurements
        pT = 1./pV;
        
        % We should be surprised if precision is decreasing
        % relative to the prior. But it might happen somtimes? unclear
        %assert(all(pT>=Ta));
        
        % we need to remove the specific contributions from each marginal
        % variable.
        sT  = pT - dP;
        sTM = pT.*pM - dPM;
        sM  = sTM./sT;
        sV  = 1./sT;
        
        % Precision should usually increase, although it <might> go down
        % in cases where multiple, conflicting pieces of information are
        % added?
        %assert(all(sT>=Ta));
        
        if rms(M-sM)/std(M)<tol && rms(sqrt(V)-sqrt(sV))/std(sqrt(V))<tol,
            % done
            %fprintf('early convergence at iteration %f\n',iteration);
            break;
        end
        
        % In theory we should now have a residual marginal that incoporates
        % information from the prior, and all other updates excluding the 
        % local update. We should be able to iterate at this point.
        V = sV;
        M = sM;
        T = sT;
    end
    
    % Iterations completed, we perform one last integration
    [pxy,M,pV,pP,dM,dP,dPM] = integrateMarginals(M,T,x,logpyx);
    
