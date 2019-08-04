function op=gaussian1DblurOperator(n,sigma,truncate)
    %{
    Returns a 1D Gaussan convolution operator of size n. The result is
    an n x n matrix that
    applies a Gaussian blur in 1D with standard deviation sigma
    when multiplied by a column vector.

    Because Gaussian is convolution is separable, 2D convolutions can
    be performed by first blurring in the x and then the y direction.
    If the blur radius (sigma) is small, and if sparse arrays are
    used, this can be of equal or greater efficiency than using 
    Fourier transforms or naive multiplication.
    
    Note that is sigma is large (i.e. blurring "a lot") then this 
    operator will be close to singular and the inverse operator
    may not exist.
    
    Parameters
    ----------
    n : int
        Should be greater than zero. 
        Length of 1D domain for which to generate the blur operator
    sigma : float
        Standard deviation of Gaussian blur operator, in units of
        pixels (samples)
    truncate : float
        Should be positive and small. Values smaller than this 
        will be rounded down to zero to make the operator sparse.
        
    Returns
    -------
    op : `matrix`
        n x n matrix which, when multiplied by a column vector, will
        apply a Gaussian blur in 1D with standard deviation sigma.

    Example
    -------
    ::
    
        % Test code for performing a 2D separable convolution
        % Size of 1D spatial domain
        n = 100; 
        % Standard deviation of blur
        sigma = 5; 
        % Construct operator
        blur = gaussian1DblurOperator(n,sigma); 

        % randomly draw some quantized (0,1) data in [0,1]²
        y = randn(n,n)<-2; 

        % Illustrate 2D blur via separable convolution
        subplot(221);
        imshow(y);
        subplot(222);
        imshow(10*blur*y);
        subplot(223);
        imshow(10*y*blur');
        subplot(224);
        imshow(10*blur*y*blur');

    %}

    % Round values smaller than `truncate` down to zero to create a 
    % sparse operator.
    if nargin<3, truncate = 1e-4; end
    
    x   = linspace(0,n-1,n); % 1D domain
    tau = 1.0/sigma^2;       % precision
    k   = exp(-tau*x.^2);    % compute (un-normalized) 1D kernel
    tp  = toeplitz(k,k);     % convert to an operator from n -> n
    op = bsxfun(@rdivide, tp', sum(tp,2))';
    op = sparsifyUnitary(op,truncate);
end
    

