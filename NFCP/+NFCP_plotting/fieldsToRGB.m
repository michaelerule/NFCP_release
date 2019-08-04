function RGBshow = fieldsToRGB(M,n,cscale,upscale,softmax)
    %{ 
    Render fields to an RGB image with upsampling
    
    Parameters
    ----------
    M: vector
        Means for concentrations of species M, concatenated
    n : scalar positive integer, typically between 5 and 50
        n x n simulation grid size.
    cscale : 3 x 1 array
        Three multipliers to convert the Q, A, R, field values to blue,
        red, and green color values, respectively
    upscale : positive integer
        Number of times to upsample image before displaying. If 1, it 
        will display at the resolution of the basis projection.
    softmax : bool, optional
        If tue, normalize the color-scheme to the peak value in each
        channel
    %}
    
    Q = M(1:n*n);
    A = M(n*n+1:n*n*2);
    R = M(n*n*2+1:n*n*3);
    
    % Normalize to 1 (locally)
    Total = Q+A+R;
    Q = Q./Total;
    A = A./Total;
    R = R./Total;
    
    b = reshape(cscale(1).*Q,n,n);
    r = reshape(cscale(2).*A,n,n);
    g = reshape(cscale(3).*R,n,n);
    RGBshow = imresize(cat(3,r,g,b),upscale);
    
    if nargin<5
        softmax = false;
    end
    if softmax,
        RGBshow = exp(5.*RGBshow);
        Total = sum(RGBshow,3);
        RGBshow = RGBshow./cat(3,Total,Total,Total);
    end
    
