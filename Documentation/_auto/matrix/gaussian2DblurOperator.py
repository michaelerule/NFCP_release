#!/usr/bin/env python
# -*- coding: UTF-8 -*-
def gaussian2DblurOperator(n,sigma,truncate):
    r'''

    Returns a 2D Gaussan blur operator for a n x n sized domain
    Constructed as a tensor product of two 1d blurs of size n.

    See gaussian1DblurOperator for an example of how to perform a
    2D Gaussian blur that is typically faster than constructing and
    applying the fully 2D blur operator (as done here).

    Note that is sigma is large (i.e. blurring "a lot") then this
    operator will be close to singular and the inverse operator
    may not exist.

    Although directly and constructing a full 2D Gaussian blur
    operator is more expensive than using the 1D separable approach,
    it is easier to define operations on the full spatiotemporal
    covariance matrix in terms of an explicit construction of the 2D
    blur operator.

    Parameters
    ----------
    n : int
        Should be greater than zero.
        Length of 2D domain for which to generate the blur operator
    sigma : float
        Standard deviation of Gaussian blur operator, in units of
        pixels (samples)
    truncate : float
        Should be positive and small. Values smaller than this
        will be rounded down to zero to make the operator sparse.

    Returns
    -------
    op : n² x n² matlab array
        n² x n² matrix. 2D data should be unraveled in row-major order
        into a n² length column vector. Multiplying op by this vector
        will apply a 2D Gaussian blur with standard deviation sigma.

    Example
    -------
    ::

        % Test code for performing a 2D separable convolution
        % Size of spatial domain
        n = 15;
        % Standard deviation of blur
        sigma = 1;
        % Construct operator
        blur = gaussian2DblurOperator(n,sigma);

        % randomly draw some quantized (0,1) data in [0,1]²
        y = randn(n,n)<-2;

        % Illustrate 2D blur
        subplot(121);
        imagesc(y);
        subplot(122);
        imagesc(reshape(blur*reshape(y,n*n,1),n,n));

    '''
    pass#SKIPME
    '''#STARTCODE

    % Round values smaller than `truncate` down to zero to create a
    % sparse operator.
    if nargin<3, truncate = 1e-4; end

    x   = linspace(0,n-1,n); % 1D domain
    tau = 1.0/sigma^2;       % precision
    k   = exp(-tau*x.^2);    % compute (un-normalized) 1D kernel
    tp  = toeplitz(k,k);     % convert to an operator from n -> n
    d2  = kron(tp,tp);       % take the tensor product to get 2D operator
    op  = bsxfun(@rdivide, d2', sum(d2,2))';
    %op  = sparsifyUnitary(op,truncate);
end



    '''#STOPCODE
