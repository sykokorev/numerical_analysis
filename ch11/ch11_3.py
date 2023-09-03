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


from integration.integrate import trapezoid, rectangle


def func(x):
    return 2 + 2 * x + x ** 2 + np.sin(2 * np.pi * x) + np.cos(2 * np.pi * x / 0.5)


def funcint(x):
    return 2 * x + x ** 2 + (1 / 3) * x ** 3 - np.cos(2 * np.pi * x) / (2 * np.pi) + np.sin(2 * np.pi * x / 0.5) / (45 * np.pi)


def exp(x):
    return np.exp(x)


if __name__ == "__main__":

    a = 0
    b = 1.5
    n = 8
    h = (b - a) / n

    xdata = np.linspace(a, b, n)
    ydata = func(xdata)
    x = np.linspace(a, b, 100)
    y = func(x)

    It = trapezoid(func, xdata)
    I0 = funcint(b) - funcint(a)
    print(f'Exact integral {I0 = }')
    print(f'Trapezoidal intgral with step {h = }, {It = }')
    print(f'Error {abs(I0 - It)}')

    fig, ax = plt.subplots()
    ax.plot(x, y, marker='', color='k', label='func')
    ax.scatter(xdata, ydata, marker='o', color='k')

    ax.grid(True)
    fig.legend()
    plt.show()

    print('Exponential function')
    a = 0.0
    b = 2.0
    n = 72
    h = (b - a) / (n - 1)
    x = np.linspace(a, b, n)
    It = trapezoid(exp, x)
    Ir = rectangle(exp, x, method='mid')
    I0 = exp(b) - exp(a)
    E = exp(b) * (b - a) * (h ** 2) / 12

    print(f'Integration by trapezoid method: {It = }')
    print(f'Integration by rectangle mid point: {Ir = }')
    print(f'Exact itegration: {I0 = }')
    print(f'The upper bound error: {E = }. The trapezoid method error: {abs(I0 - It)}')

