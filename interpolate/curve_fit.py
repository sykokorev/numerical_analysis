import os
import sys
import math


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(DIR)


from la.linalg import *


def simple_linear_model(a, b, x):
    return a * x + b


def linear_model(x: float, f: list, c: list) -> float:
    if len(f) != len(c): raise ValueError('Length of the list of funciton has to be equal length of the list of coefficients.')
    return sum([ci * fi(x) for ci, fi in zip(c, f)])


def exp_model(a, b, x):
    return b * math.e ** (a * x)


def power_model(a, b, x):
    return b * x ** a


def jacobian_row(f, x, args, del_args=1e-8):

    n = len(args)
    j = zeros(n)

    for i in range(n):

        tmp = [args[j] if i != j else args[j] + del_args for j in range(n)]
        j[i] = f(x, *args) - f(x, *tmp)
    return j


def jacobian(f, x, args, del_args=1e-8):

    if not isinstance(f, list):
        if not hasattr(f, '__call__'):
            raise ValueError('First parameter must be function or list of functions')
        else:
            f = [f]
    if not any([hasattr(fi, '__call__') for fi in f]):
        raise ValueError('First parameter must be function or list of functions')
    
    m = len(f)
    n = len(args)
    j = [[0.0 for i in range(n)] for j in range(m)]
    for i in range(m):
        j[i] = jacobian_row(f[i], x, args, del_args)
    
    return j[0] if m == 1 else j


def simple_linear(x: list, y: list) -> tuple:

    '''
        Linear regression model y = ax + b
        Returns tuple of coefficients (a, b) and correlation coefficient R^2
    '''

    if len(x) != len(y):
        raise ValueError('x and y data must have equal length.')

    n = len(x)
    x_sum = sum([xi for xi in x])
    y_sum = sum([yi for yi in y])
    xsq_sum = sum([xi ** 2 for xi in x])
    xy_sum = sum([xi * yi for xi, yi in zip(x, y)])
    y_ave = y_sum / n

    a = (n * xy_sum - x_sum * y_sum) / (n * xsq_sum - x_sum ** 2)
    b = (y_sum - a * x_sum) / n
    
    s1 = sum([(yi - simple_linear_model(a, b, xi)) ** 2 for xi, yi in zip(x, y)])
    s2 = sum([(yi - y_ave) ** 2 for yi in y])

    r_sq = 1 - s1 / s2

    return ((a, b), r_sq)


def ext_linear(f: list, x: list, y: list) -> tuple:

    '''
        --------------------------------------------------------------------------------------
        Extention of a linear regression model
        Given a set of data points and a model y(x) = a1 * f1(x) + a2 * f2(x) + ... am * fm(x)
        Find a1, a2, ..., am so that the model best fits to the set of data points
        --------------------------------------------------------------------------------------
        f is a list of m functions
        x and y are the lists of m data points
        return tuple of a coefficients a1, a2, ..., am and a correlation coefficient
    '''

    n = len(f)
    if len(x) != len(y): raise ValueError('x and y data points have to be equal length')

    y_ave = sum(y) / len(y)

    fm = zeros(n, axis=2)
    ym = zeros(n, axis=1)

    for i, fi in enumerate(f):
        fm[i] = [sum([fi(xi) * fj(xi) for xi in x]) for fj in f]
        ym[i] = sum([yi * fi(xi) for xi, yi in zip(x, y)])

    coef = solveLUP(fm, ym)
    s1 = sum([(yi - linear_model(xi, f, coef)) ** 2 for xi, yi in zip(x, y)])
    s2 = sum([(yi - y_ave) ** 2 for yi in y])

    r_sq = 1 - s1 / s2

    return (coef, r_sq)


def exp(x: list, y: list) -> tuple:

    v = [math.log(yi) for yi in y]
    u = [xi for xi in x]

    (c1, c2), r =  simple_linear(u, v)
    b, a = math.e ** c2, c1
    s1 = sum([(yi - exp_model(xi, b, a)) ** 2 for xi, yi in zip(x, y)])
    s2 = sum([yi ** 2 for yi in y])
    r = 1 - s1 / s2

    return (a, b), r


def power(x, y):

    v = [math.log(yi) for yi in y]
    u = [math.log(xi) for xi in x]

    (c1, c2), r = simple_linear(u, v)
    a, b = c1, math.e ** c2

    s1 = sum([(yi - power_model(a, b, xi)) ** 2 for xi, yi in zip(x, y)])
    s2 = sum([yi ** 2 for yi in y])
    r = 1 - s1 / s2

    return (a, b), r
