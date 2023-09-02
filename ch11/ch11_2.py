import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt


from scipy.optimize import curve_fit


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(DIR)


from integration.integrate import rectangle, polyint


def func(x):
    return x ** 2


def exp(x):
    return np.exp(x)


if __name__ == "__main__":

    h = 4   
    a = 0.0
    b = 1.0
    x = np.linspace(a, b, h)
    y = func(x)
    I0 = polyint(b, 0, 0, 1) - polyint(a, 0, 0, 1)
    I1 = rectangle(func, x, method='mid')
    I2 = rectangle(func, x, method='right')
    I3 = rectangle(func, x, method='left')

    print(I0, I1, I2, I3)

    a = 0.0
    b = 2.0
    h = 0.004
    h = int((b - a) / h)
    x = np.linspace(a, b, h)
    y = exp(x)
    
    I0 = exp(b) - exp(a)
    I1 = rectangle(exp, x, method='left')
    I2 = rectangle(exp, x, method='mid')
    I3 = rectangle(exp, x, method='right')

    print(f'Explicit {I0=}')
    print(f'Left sum(f(xi-1)), {I1=}')
    print(f'Mid sum(f((xi + xi-1) / 2)), {I2=}')
    print(f'Right sum(f(xi)), {I3=}')

    h = (b - a) / h
    print('Error')
    e1 = exp(b) * (b - a) * h / 2
    e2 = exp(b) * (b - a) * (h ** 2) / 24
    e3 = exp(b) * (b - a) * h / 2
    print(f'Upper bound left error: {e1}, error: {abs(I0 - I1)}')
    print(f'Upper bound mid error: {e2}, error: {abs(I0 - I2)}')
    print(f'Upper bound right error: {e3}, error: {abs(I0 - I3)}')
