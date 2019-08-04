function chC=robustChol(C)
    %{
    Cholesky factor, but handle non-PSD matrices
    
    Parameters
    ----------
    x : matrix
        Positive definite square matrix
    
    Returns
    -------
    matrix
        The inverse of x
    %}
    try
        chC = chol(C);
    catch
        try
            chC = chol(0.5*(C+C')); % x = chol(x)'*chol(x)
        catch
            chC = chol(nearestSPD(C));
        end
    end
