function x=cinv(x)
    %{
    Fast inverse of positive matrices using Cholesky decomposition.
    This routine will fail (when trying to compute the Cholesky
    decomposition) if the matrix is not positive definite.
    
    Parameters
    ----------
    x : matrix
        Positive definite square matrix
    
    Returns
    -------
    matrix
        The inverse of x
    %}
    ch = chol(0.5*(x+x')); % x = chol(x)'*chol(x)
    ch = ch\eye(size(ch,1));
    x  = ch*ch';
