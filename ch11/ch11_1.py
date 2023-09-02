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


from interpolate.optimize import lm, polyfit
from diff.diff import derivative as der
from integration.integrate import polyint


def func(x, a, b, c):
    return a - x ** 2 + b * np.cos(2 * math.pi * x / c)


def funcint(x):
    return 2 * x - (x ** 3) / 3 + 0.1 * (0.7 / (2 * math.pi)) * np.sin(2 * math.pi * x / 0.7)


def runge(x):
    return 1 / (1 + 25 * x ** 2)


def poly(x, *args):
    return sum([c * x ** i for i, c in enumerate(args)])


if __name__ == "__main__":

    args = [2, 0.1, 0.7]
    a, b = -1.0, 1.0 # Integration limits
    h = 21 # Number of points divided the function func

    x = np.linspace(a, b, 100)
    y = func(x, *args)

    xdata = np.linspace(a, b, h)
    popt = polyfit(x, y, h)
    print(popt)
    ydata = poly(xdata, *popt)

    Ia = funcint(a)
    Ib = funcint(b)
    print(f'Integral: {Ib - Ia}')

    Fb = polyint(b, *popt)
    Fa = polyint(a, *popt)
    print(f'Approximated polynomial order {h - 1} integral = {Fb - Fa}')
    print(f'Absolute error: {abs((Ib - Ia) - (Fb - Fa))}')

    fig, ax = plt.subplots()
    ax.plot(x, y, marker='', color='k', label='func')
    ax.scatter(xdata, ydata, marker='o', color='r', label='polynomial')
    ax.plot(xdata, ydata, marker='', color='g')

    ax.grid(True)
    fig.legend()
    
    plt.show()
