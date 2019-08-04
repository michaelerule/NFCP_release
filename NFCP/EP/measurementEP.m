function [loglikelihood,M2,C2,Ma2,Pa2,rMa,rPa,rCa] = measurementEP(model,M,C,xypoints,srMa,srPa,srCa)
    %{
    Solve non-conjugtae Gaussian-Poisson update using expectation
    propagation. 
    
    Parameters
    ----------
    model : the model structure

    %}
    m    = model.nn;   % No. basis function sper state
    Ir   = eye(m).*model.reg_inverse;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract active state
    assertPSD(C,model.safety);
    Ma   = M(m+1:2*m,1);
    Ca   = C(m+1:2*m,m+1:2*m);
    assertPSD(Ca,model.safety);
    chCa = chol(Ca); % x = chol(x)'*chol(x)
    chPa = (chCa\eye(m))';
    Pa   = chPa'*chPa;
    assertPSD(Pa,model.safety);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Observations
    y    = binCounts(model,xypoints);   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EP to get directly observed subspace
    [pMa,pPa] = expectationPropagationUpdate(model, Ma, Pa, Ca, y);
    assertPSD(pPa,model.safety);
    
    % EP to get bias from domain truncation (removed negative values)
    [BMa,BPa] = expectationPropagationUpdate(model, Ma, Pa, Ca, y.*NaN);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % log-likelihood
    loglikelihood = expectedLogLikelihood(model,Ma,Pa,pMa,pPa,y);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Kalman update to get joint
    [M2,C2,rPa,rCa,rMa] = KalmanPropagate(M,C,Pa,Ma,pPa+Ir,pMa);
    assertPSD(C2,model.safety);
    
    Ma2 = pMa;
    Pa2 = pPa;
