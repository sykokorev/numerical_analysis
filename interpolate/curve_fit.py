import os
import sys


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
