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


def lm(f, xdata, ydata, p0=None, jac=None, lam0 = 10, est = 1e-6, iter=0):

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
                lam0: float
                    Initial damping factor
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

    n = len(xdata)
    nargs = len(signature(f).parameters) - 1
    popt = [1.0 for i in range(nargs)] if not p0 else p0
    delarg = 1e-6
    J = [[0.0 for i in range(nargs)] for j in range(n)]
    fx = [f(xi, *popt) for xi in xdata]
    
    if not jac:
        for i in range(n):
            for j in range(nargs):
                args1 = [popt[k] + delarg * popt[k] if k == j else popt[k] for k in range(nargs)]
                args2 = [popt[k] - delarg * popt[k] if k == j else popt[k] for k in range(nargs)]
                J[i][j] = (f(xdata[i], *args1) - f(xdata[i], *args2)) / (2 * delarg * popt[j])
    else:
        for i in range(n):
            J[i] = jac(xdata[i], *popt)

    JT = transpose(J)
    H = mat_mul(JT, J)
    H_diag = [[H[i][j] if i == j else 0.0 for i in range(nargs)] for j in range(nargs)]
    H_lam = mat_add(H, mat_mul(lam0, H_diag))
    rhs = mat_mul(JT, [yi - fxi for yi, fxi in zip(ydata, fx)])
    dp = mat_mul(rhs, inverse(H_lam))

    popt_new = [p + dpi for p, dpi in zip(popt, dp)]
    fx_new = [f(xi, *popt_new) for xi in xdata]

    rmse = (sum([(yi - fi) ** 2 for yi, fi in zip(ydata, fx)]) / n)
    rmse_new = (sum([(yi - fi) ** 2 for yi, fi in zip(ydata, fx_new)]) / n)

    es = abs(rmse - rmse_new) / rmse

    if es < est:
        return popt_new
    elif rmse_new > rmse:
        return lm(f, xdata, ydata, popt, jac, lam0 * 10, est, iter)
    else:
        return lm(f, xdata, ydata, popt_new, jac, lam0 / 10, est, iter+1)
