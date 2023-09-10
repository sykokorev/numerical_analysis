import os
import sys
import numpy as np
from math import floor


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(DIR)


from la.linalg import zeros


def factor(n):

    if n == 0:
        return 1
    elif n == 1:
        return 1
    else:
        return n * factor(n - 1)
    

def binomial(n, k):
    if n < k:
        return 0
    return factor(n)  / (factor(k) * factor(n - k))


def legendre_poly_coef(n: int) -> list:
    '''
        Compute coefficients of Legendre polynomials degree n
        -----------------------------------------------------
            Paramters: n: int
                Polynomial degree
            Returns: list
                returns list of coefficients and degree of undependet veriable
                of the Legendre polynomial
                [(a0, k0), (a1, k1), (a2, k2), (a3, k3) ...]
                L(x) = a0 * x ** k + a1 * x ** k1 + a2 * x ** k2 + ...
    '''
    const = 1 / 2 ** n
    return [(
        const * (-1) ** k * binomial(n, k) * binomial(2 * n - 2 * k, n), 
        n - 2 * k
        ) for k in range(0, floor(n / 2) + 1)]


def legendre_poly(x, *args):

    return sum([a[0] * x ** a[1] for a in args])


def polyint(x, *coeff):
    return sum([(1 / (i + 1)) * c * x ** (i + 1) for i, c in enumerate(coeff)])


def rectangle(f: callable, a: int | float, b: int | float, k: int = 2, args = [], method = 'mid'):

    '''
        Compute integral by using rectangle method
        ------------------------------------------
            f: callable: Function that need be integrated
            a, b: int | float: Integration limits
            k: int (default 2): Number of intervals
            args: list (default []): Extra arguments for function f(x, *args)
            method: str (default 'mid'): Method that will be used ('mid', 'right', 'left')
    '''

    I = 0
    x = np.linspace(a, b, k + 1)
    h = abs(x[0] - x[1])

    if method == 'mid':
        I = h * sum([f((x[i - 1] + x[i]) / 2, *args) for i in range(1, k + 1)])
    elif method == 'right':
        I = h * sum([f(x[i], *args) for i in range(1, k + 1)])
    elif method == 'left':
        I = h * sum([f(x[i - 1], *args) for i in range(1, k + 1)])

    return I


def trapezoid(f: callable, a: int | float, b: int | float, k: int = 2, args: list = []):

    '''
        Compute intgral by using trapezoid rule
        ---------------------------------------
            f: callable: function which intgral will be computed
            a, b: int | float: intagration limits
            k: int (default 2): number of parts on which interval will be divided
            args: list (default []): other function arguments f(x, *args)
        ---------------------------------------
            return value of integration and absolute error equal |E| = abs(Ik - Ik-1)
        
    '''

    h = (b - a) / k
    I = f(a, *args) + f(b, *args)

    for i in range(1, k):
        xi = a + i * h
        I = I + 2 * f(xi, *args)
    
    I = I * h / 2
    return I


def simpson(f: callable, a: int | float, b: int | float, k: int = 2, args = [], method=1):

    '''
        Compute the integral of the function f on the interval [a, b]
        -------------------------------------------------------------
            f: callable: function which integral will be computed
            a, b: int | float: Integration limits
            k: int (default 2): Number of subdivision
            args: list (default []): other function arguments f(x, *args)
            method: int (default 1): method with which integral will be computed
                    1 - Simpson's 1/3 rule
                    2 - Simpson's 3/8 rule
        -------------------------------------------------------------
            return value of integration and absolute error |E| = abs(Ik - Ik-1)
    '''

    x = np.linspace(a, b, k + 1)
    I = 0.0
    if method == 1:
        h = (x[1] - x[0]) / 2
        for i in range(1, k + 1):
            if i == k:
                Ik1 = I
            xm = x[i - 1] + h
            I += (f(x[i - 1], *args) + 4 * f(xm, *args) + f(x[i], *args)) * (h / 3)
    if method == 2:
        h = (x[1] - x[0]) / 3
        for i in range(1, k + 1):
            if i == k:
                Ik1 = I
            xli = x[i - 1] + h
            xri = x[i - 1] + 2 * h
            I += (f(x[i - 1], *args) + 3 * f(xli, *args) + 3 * f(xri, *args) + f(x[i])) * (3 * h / 8)

    e = abs(I - Ik1)

    return I, e


def romberg(f: callable, a: float, b: float, args: list = (), tol=1.48e-08, rtol=1.48e-08) -> float:

    '''
        Romberg integration of a callable function. Returns the integral of a function
        ------------------------------------------------------------------------------
        Parameters:     f: callable: 
                            function to be integrated
                        a, b: float: 
                            Lower and upper integration limits
        Returns:        results: float
                            Result of the integration

        Other           args: tuple, optional
                            Extra argument to pass to function. Default is to pass no extra aguments
                        tol, rtol: float
                            The desired absolute and relative tolerances. Default are 1.48e-08
    '''

    # Initial setup maximum number of iterations 800
    max_iter = 800
    # First Romberg's matrix is 2x2 size
    r = zeros(2, 2)

    for iter in range(max_iter):

        # Increase matrix size if current iteration is not first. Old values are retained
        if iter:
            for k in range(iter + 2):
                if k <= iter:
                    r[k].append(0.0)
                else:
                    r.append([0.0 for i in range(iter + 2)])
        
        for i in range(iter + 2):
            
            # Fill matrix at first step. If current iteration is not first, then add additional values
            if iter == 0:
                for j in range(iter, iter - i + 2):
                    if i == 0:
                        r[i][j] = trapezoid(f, a, b, 2 ** j)
                    elif j < iter + 1:
                        const = 2 ** (2 * i)
                        r[i][j] = (const * r[i - 1][j + 1] - r[i - 1][j]) / (const - 1)
            else:
                for j in range(iter + 1 - i, iter + 2):
                    if i == 0:
                        r[i][j] = trapezoid(f, a, b, 2 ** j)
                    elif j == (iter + 1 - i):
                        const = 2 ** (2 * i)
                        r[i][j] = (const * r[i - 1][j + 1] - r[i - 1][j]) / (const - 1)
        
        es = abs(r[-1][0] - r[-2][0])
        er = es / r[-1][0]

        if es <= tol or er <= rtol:
            return r[-1][0]
    
    return None


def gauss_quad(f: callable, n: int = 2) -> float:
    '''
        Gaussian quadrature intgration of callable function. 
        Return the integral of a function f
        ----------------------------------------------------
        Parameters:     f: callable
                            function to be integrated
                        n: number of points to approximate
        Returs:         results: float
                            result of integration
    '''

    pass
