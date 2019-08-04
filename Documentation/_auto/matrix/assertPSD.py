#!/usr/bin/env python
# -*- coding: UTF-8 -*-
def assertPSD(H,safety,tolerance):
    r'''

    If a precision matrix is near-singular, we cannot find a gradient.
    This checks the reciprocal condition number of the precision
    or Hessian matrix. If it is smaller than the tolerance, an error
    is raised.

    Parameters
    ----------
    H : matrix
        Precision or Hessian matrix from Newton-Raphson update
        for Poisson likelihood
    safety : int; default 2
        Safety flag
          0 = Do nothing (skip all checks)
          1 = Repair errors when possible (this function does not repair)
          2 = Test and warn on unexpected numerical situations
          3 = Test and Hard-fail on unexpected numerical situations
    tolerance : scalar; optional, default `eps`
        Minimum reciprocal condition number permitted. If rcond(H)
        is smaller than this, an error will be thrown. Defaults to
        `eps`.

    '''
    pass#SKIPME
    '''#STARTCODE
    if nargin<2,
        safety=2;
    end;
    if nargin<3,
        tolerance=eps;
    end;

    % Skip checks for low safety levels
    if safety>=2,
        assertMatrixCondition(H,safety,tolerance);
        try
            ch = chol(H);
        catch
            H = nearestSPD(H);
            if safety==3,
                error('Matrix is not positive semidefinite');
            else
                warning('Matrix is not positive semidefinite');
            end
        end
    end

    '''#STOPCODE
