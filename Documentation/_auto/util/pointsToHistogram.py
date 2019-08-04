#!/usr/bin/env python
# -*- coding: UTF-8 -*-
def pointsToHistogram(xypoints,n,interpolate,verbosity):
    r'''

    Bin point xypoints into a 2D histogram.
    xypoints should be a list of (x,y) points.

    If it is not in this format, this function does nothing and
    returns the original xypoints array. This is to support transparently
    running on xypoints that have already been binned through some other
    pre-processing step.

    Parameters
    ----------
    xypoints:
        Either a npoints x 2 list of (x,y) points, or an n x n
        histogram. xypoints should lie in the [0,1]x[0,1] unit square.
    n:
        Dimension of grid onto which to bin. The result will be
        an n x n 2d histogram with n² entries. If xypoints are already
        binned (i.e. the xypoints argument is a square array rather than
        a npoints x 2 list of (x,y) locations), this function will
        ensure that the xypoints are n x n
    interpolate: bool, optional (default true)
        Use interpolation rather than binning to compute the histogram
        For example, if a point lands on the corner of four histogram
        bins, it's density will be split equally between the four
        bins rather than being randomly assigned to one of them. This
        behavior is similar to interp2d. Interpolation is slower than
        histogram binning.
    verbose: int, optional, default 0
        Verbosity flag
        0 = quiet
        1 = selected messages per time-step
        2 = detailed logging, once per time-step
        3 = detailed messages, including Laplace update convergence iterations

    Returns
    -------
    hist:
        1 x n² binned histogram, unraveled in row-major order,
        and packed as a row-vector.

    Example
    -------
    ::

        % Compare binned vs interpolated histograms
        K = 15;  % number of random points
        n = 20;  % number of bins along each dimension
        xy = rand(K,2); % sample random points
        subplot(121);
        imagesc(points_to_histogram(xy,n,false));
        subplot(122);
        imagesc(points_to_histogram(xy,n,true));


    '''
    pass#SKIPME
    '''#STARTCODE

    if nargin<3
        % Set to true to use interpolation when computing histogram
        interpolate = true;
    end
    if nargin<4
        % Set to true to print detailed debugging information while
        % running.
        verbosity = 0;
    end
    if numel(xypoints)<=0,
            hist = zeros(n,n);
    elseif size(xypoints,2)~=2
        if verbosity>2,
            display('Assuming xypoints provided as histogram)');
            size(xypoints)
        end
        hist = xypoints;
    else
        if verbosity>2,
            display('Assuming xypoints provided as (x,y) points');
        end
        if (any(any(xypoints<0.0)) || any(any(xypoints>1.0))),
            error('points_to_histogram: point outside [0,1]² square');
        end
        if interpolate,
            % Use binlinear interpolation.
            % 1e-9 subtracted to address edge case of x or y = 1
            %x   = xypoints(1:end,1).*(n-1-1e-9)+1;
            %y   = xypoints(1:end,2).*(n-1-1e-9)+1;
            x   = xypoints(1:end,1).*n+0.5;
            y   = xypoints(1:end,2).*n+0.5;
            x   = min(max(x,1+1e-9),n-1e-9);
            y   = min(max(y,1+1e-9),n-1e-9);
            ix  = floor(x);
            iy  = floor(y);
            fx  = x-ix;
            fy  = y-iy;
            p22 = fx.*fy;
            p21 = fx.*(1-fy);
            p12 = fy.*(1-fx);
            p11 = (1-fx).*(1-fy);
            npoints = size(xypoints,1);
            hist    = zeros(n,n);
            for j=1:npoints,
                jx = ix(j);
                jy = iy(j);
                if n>1,
                    hist(jx  ,jy  ) = hist(jx  ,jy  ) + p11(j);
                    hist(jx  ,jy+1) = hist(jx  ,jy+1) + p12(j);
                    hist(jx+1,jy  ) = hist(jx+1,jy  ) + p21(j);
                    hist(jx+1,jy+1) = hist(jx+1,jy+1) + p22(j);
                else
                    % Handle zero-dimensional (nonspatial) as special case
                    hist = hist + numel(x);
                end
            end
            hist = hist';
        else
            % Use nearest-neighbor binning
            edges{1} = 0:1/n:1;
            edges{2} = 0:1/n:1;
            hist = hist3(xypoints,'Edges',edges);
            hist = hist(1:n,1:n);
            hist = hist';
        end
    end



    '''#STOPCODE
