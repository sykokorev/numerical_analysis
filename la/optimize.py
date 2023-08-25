import os
import sys
import numpy as np


from inspect import signature


IMP = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(IMP)


from la.linear_algebra import *
from la.linalg import solveLUP as lup


def lm(f, xdata, ydata, p0=None, jac=None, est=10e-6):

    '''

    Parameters: f: callable
                    The model function, f(x, ...). It must take an independent
                    variable as the first argument and the parameters to fit as
                    separate remainin arguments.
                xdata: array_like
                    The independent variables a length m
                ydata: array_like
                    The dependent variable, a length m
                p0: array_like, optional
                    Initial guess for the parameters (length n). If None (default),
                    then the initial parameters will all be set to 1
                jac: callable, optional
                    Function with signature jac(x, ...) which computes the 
                    Jacobian matrix of the model function with respect to 
                    parameters. If None (default), the Jacobian will be estimated
                    numerically
                est: float | array_like
                    Approxiamtion error. If float then the estimation error apply to all
                    parameters of the model
    Returns: potp: array
                    Optimal values for the parameters so that the sum of the squared
                    residuals of f(xdata, *popt) - ydata is minimized
    
    '''

    if not any([isinstance(xi, float | int) for xi in xdata]):
        raise ValueError('Incompatable data for xdata values')
    if not any([isinstance(yi, float | int) for yi in ydata]):
        raise ValueError('Incompatable data for ydata values')
    if not hasattr(f, '__call__'):
        raise ValueError('First argument must be callable')

    n = len(xdata)
    nargs = len(signature(f).parameters) - 1
    popt = [1.0 for i in range(nargs)] if not p0 else p0
    delarg = 10e-6
    J = [[0.0 for i in range(nargs)] for j in range(n)]
    fx = [f(xi, *popt) for xi in xdata]

    # Compute Jacobian matrix
    if not jac:
        for i in range(n):
            for j in range(nargs):
                args1 = [popt[k] + delarg if k == j else popt[k] for k in range(nargs)]
                args2 = [popt[k] - delarg if k == j else popt[k] for k in range(nargs)]
                J[i][j] = (f(xdata[i], *args1) - f(xdata[i], *args2)) / (2 * delarg)
    else:
        for i in range(n):
            J[i] = jac(xdata[i], *popt)

    # Compute Hessian and gradient of the function
    JT = transpose(J)
    H = mat_mul(JT, J)
    g = mat_mul(JT, fx)

    # Compute step
    step = lup(H, mat_mul(JT, [yi - f(xi, *popt) for xi, yi in zip(xdata, ydata)]))
    popt_new = [p + pn for p, pn in zip(popt, step)]

    # Check estimation error
    if all([abs(p - pn) <= est for p, pn in zip(popt, popt_new)]):
        return popt_new, [abs(p - pn) for p, pn in zip(popt, popt_new)]
    else:
        return lm(f, xdata, ydata, popt_new, jac, est)
