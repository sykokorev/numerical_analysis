import os
import sys
import numpy as np


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__), 
    os.pardir
))
sys.path.append(DIR)

from la.linalg import newton_raphson2



def euler(func, t_span, y0, h, args=()) -> tuple:

    '''
        Solves ODE IVP by using (explicit) Euler method
        -----------------------------------------------
        Paramters:      func: callable
                            Right-hand side of the equation. 
                            The time derivative of the state y at time t
                            func(x, t, *args)
                        t_span: 2-members sequence
                            Interval of integration (t0, tf)
                        y0: float
                            Initial state
                        steps: int
                            Number of integration steps
                        args: tuple, optional
                            Extra arguments to pass to the function
    '''

    if hasattr(y0, '__iter__'):
        ysol = [np.asanyarray(y0)]
        def fun(t, x, func=func, *args):
            return np.asarray(func(t, x, *args))
    else:
        ysol = [y0]
        def fun(t, x, func=func, *args):
            return func(t, x, *args)

    tsol = [t_span[0]]
    n = (int)((t_span[1] - t_span[0]) / h)
    
    for i in range(n):
        if tsol[i] + h <= t_span[1]:
            ysol.append(ysol[i] + h * fun(tsol[i], ysol[i], func, *args))
            tsol.append(tsol[i] + h)

    if abs(tsol[-1] - t_span[1]) > 1e-3:
        h = t_span[1] - tsol[-1]
        ysol.append(ysol[-1] + h * fun(tsol[-1], ysol[-1], func, *args))
        tsol.append(t_span[1])

    return tsol, ysol


def heun(func, t_span, y0, h, args=()) -> tuple:

    '''
        Solves ODE IVP by using Heun's method
        -----------------------------------------------
        Paramters:      func: callable
                            Right-hand side of the equation. 
                            The time derivative of the state y at time t
                            func(x, t, *args)
                        t_span: 2-members sequence
                            Interval of integration (t0, tf)
                        y0: float
                            Initial state
                        steps: int
                            Number of integration steps
                        args: tuple, optional
                            Extra arguments to pass to the function
    '''

    if hasattr(y0, '__iter__'):
        ysol = [np.asanyarray(y0)]
        def fun(t, x, func=func, *args):
            return np.asarray(func(t, x, *args))
    else:
        ysol = [y0]
        def fun(t, x, func=func, *args):
            return func(t, x, *args)

    tsol = [t_span[0]]
    n = int((t_span[1] - t_span[0]) / h)

    for i in range(n):

        if tsol[i] + h <= t_span[1]:
            y_ = ysol[i] + h * fun(tsol[i], ysol[i], func, *args)
            ysol.append(ysol[i] + h * (fun(tsol[i] + h, y_, func, *args) + fun(tsol[i], ysol[i], func,*args)) / 2)
            tsol.append(tsol[i] + h)

    if abs(tsol[-1] - t_span[1]) > 1e-3:
        h = t_span[1] - tsol[-1]
        y_ = ysol[-1] + h * fun(tsol[-1], ysol[-1], func, *args)
        ysol.append(ysol[i] + h * (fun(tsol[-1] + h, y_, func, *args) + fun(tsol[-1], ysol[-1], func, *args)) / 2)
        tsol.append(tsol[-1] + h)
    
    return tsol, ysol


def mid_point(func, t_span, y0, h, args=()) -> tuple:

    '''
        Solves ODE IVP by using mid point method
        -----------------------------------------------
        Paramters:      func: callable
                            Right-hand side of the equation. 
                            The time derivative of the state y at time t
                            func(x, t, *args)
                        t_span: 2-members sequence
                            Interval of integration (t0, tf)
                        y0: float
                            Initial state
                        steps: int
                            Number of integration steps
                        args: tuple, optional
                            Extra arguments to pass to the function
    '''

    if hasattr(y0, '__iter__'):
        ysol = [np.asanyarray(y0)]
        def fun(t, x, func=func, *args):
            return np.asarray(func(t, x, *args))
    else:
        ysol = [y0]
        def fun(t, x, func=func, *args):
            return func(t, x, *args)

    tsol = [t_span[0]]
    n = (int)((t_span[1] - t_span[0]) / h)

    for i in range(n):
        if tsol[i] + h <= t_span[1]:
            y_ = ysol[i] + h * fun(tsol[i] , ysol[i], func, *args) / 2
            ysol.append(ysol[i] + h * fun(tsol[i] + h / 2, y_, func, *args))
            tsol.append(tsol[i] + h)

    if abs(tsol[-1] - t_span[1]) > 1e-3:
        h = t_span[1] - tsol[-1]
        y_ = ysol[-1] + h * fun(tsol[-1], ysol[-1], func, *args) / 2
        ysol.append(ysol[-1] + h * fun(tsol[-1] + h / 2, y_, func, *args))
        tsol.append(tsol[-1] + h)

    return tsol, ysol


def trapz(func, t_span, y0, h, args=()):

    if hasattr(y0, '__iter__'):
        ysol = [np.asanyarray(y0)]
        def fun(t, x, func=func, *args):
            return np.asarray(func(t, x, *args))
    else:
        ysol = [y0]
        def fun(t, x, func=func, *args):
            return func(t, x, *args)

    tsol = [t_span[0]]
    n = (int)((t_span[1] - t_span[0]) / h)

    for i in range(n):
        if tsol[i] + h <= t_span[1]:
            f = fun(tsol[i] + h, ysol[i], func, *args) + fun(tsol[i], ysol[i], func, *args)
            ysol.append(ysol[i] + 0.5 * h * f / (1 - 0.5 * h * fun(tsol[i], ysol[i] + 1e-6, func, *args)))
            tsol.append(tsol[i] + h)
    
    if abs(tsol[-1] - t_span[1]) > 1e-3:
            h = t_span[1] - tsol[-1]
            f = fun(tsol[-1] + h, ysol[-1], func, *args) + fun(tsol[-1], ysol[-1], func, *args)
            ysol.append(ysol[-1] + 0.5 * h * f / (1 - 0.5 * h * fun(tsol[-1], ysol[-1] + 1e-6, func, *args)))
            tsol.append(tsol[-1] + h)

    return tsol, ysol


def solve_ipv(fun, tspan, ic, method='RK45', teval=None, args=()):
    '''
        Solve an initial value problem for ODE
        --------------------------------------
        Parameters:     fun: callable
                            Right-hand side of the ODE: the time derivative
                            of the state y at time t. The calling signature
                            func(t, y), where t is a scalar and y is an 
                            ndarray with len(y) = len(ic). func must return
                            an array of the same shape as y.
                        tspan: 2-member sequence
                            Interval of integration (t0, tf)
                        ic: float or array_like, shape (n,)
                            Initial conditions.
                        intervals: integer, optional
                            Number of intervals. Default 100
                        method: string, optional
                            Integration method to use:
                                *   'RK45' (default): Explicit Runge-Kutta method of
                                    order 4
                        teval: array_like or None, optional
                            Times at which to store the computed solution, must be lie
                            within tspan. If None (default), use points selected by solver
        Rerturns:       t: ndarray, shape(n_points)
                        y: ndarray, shape(n, n_points)
                            
    '''

    yf = [ic] if not hasattr(ic, '__iter__') else ic
    n = 0 if not hasattr(ic, '__iter__') else len(ic)
    tf = tspan[0]
    tsol = [tf]
    ysol = [[y0i] for y0i in yf]
    teval = [tspan[0] + i * (tspan[1] - tspan[0]) / 100 for i in range(100)] if not hasattr(teval, '__iter__') else teval

    def f(t, y, func, *args):
        if n == 0:
            return func(t, *y, *args)
        else:
            return func(t, y, *args)
    
    for t in range(len(teval) - 1):
        h = teval[t + 1] - teval[t]
        k1 = [f(tf, yf, fun, *args)] if n == 0 else f(tf, yf, fun, *args)
        k2 = [f(tf + h / 2, [yfi + h * k1i / 2 for yfi, k1i in zip(yf, k1)], fun, *args)] if n == 0 \
            else f(tf + h / 2, [yfi + h * k1i / 2 for yfi, k1i in zip(yf, k1)], fun, *args)
        k3 = [f(tf + h / 2, [yfi + h * k2i / 2 for yfi, k2i in zip(yf, k2)], fun, *args)] if n == 0 \
            else f(tf + h / 2, [yfi + h * k2i / 2 for yfi, k2i in zip(yf, k2)], fun, *args)
        k4 = [f(tf + h, [yfi + h * k3i for yfi, k3i in zip(yf, k3)], fun, *args)] if n == 0 \
            else f(tf + h, [yfi + h * k3i for yfi, k3i in zip(yf, k3)], fun, *args)
        
        yf = [
            yfi + (h / 6) * (k1i + 2 * k2i + 2 * k3i + k4i) 
            for yfi, k1i, k2i, k3i, k4i in zip(yf, k1, k2, k3, k4)
        ]
        tf += h
        tsol.append(tf)
        if n == 0:
            ysol.append(yf)
        else:
            for i in range(n):
                ysol[i].append(yf[i])

    if n == 0:
        ysol = [yi[0] for yi in ysol]

    return tsol, ysol
