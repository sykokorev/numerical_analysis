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


from integration.integrate import trapezoid, rectangle, simpson


def func(x):
    return 2 + 2 * x + x ** 2 + np.sin(2 * np.pi * x) + np.cos(2 * np.pi * x / 0.5)


def funcint(x):
    return 2 * x + x ** 2 + (1 / 3) * x ** 3 - np.cos(2 * np.pi * x) / (2 * np.pi) + np.sin(2 * np.pi * x / 0.5) / (45 * np.pi)


def exp(x):
    return np.exp(x)


if __name__ == "__main__":

    a = 0.0
    b = 1.5
    n = 4
    h = (b - a) / (2 * (n - 1))

    x = np.linspace(a, b, 100)
    y = func(x)

    xdata = np.linspace(a, b, n)
    ydata = func(xdata)

    I0 = funcint(b) - funcint(a)
    Is1 = simpson(func, a, b, n, method=1)
    Is2 = simpson(func, a, b, n, method=2)

    print('Integral func')
    print(f'Exact integral: {I0 = }')
    print(f'Simpson\'s 1/3 rule with {h = }, {Is1[0] = }')
    print(f'Simpson\'s 3/8 rule with {h = }, {Is2[0] = }')
    print(f'Simpson\'s 1/3 error E = {abs(I0 - Is1[0])}')
    print(f'Simpson\'s 3/8 error E = {abs(I0 - Is2[0])}')

    print('Integral exp')
    a = 0.0
    b = 2.0
    n = 3
    h1 = (b - a) / (2 * (n - 1))
    h2 = (b - a) / (3 * (n - 1))
    
    x = np.linspace(a, b, 100)
    y = exp(x)

    xdata = np.linspace(a, b, n)
    ydata = exp(xdata)

    I0 = exp(b) - exp(a)
    Is1 = simpson(exp, a, b, n, method=1)
    Is2 = simpson(exp, a, b, n, method=2)
    E = exp(b) * (b - a) * (h ** 4 / 180)
    print(f'Exact integral: {I0 = }')
    print(f'Simpson\'s rule 1/3 with {h1 = }: {Is1[0] = }')
    print(f'Simpson\'s 8/3 rule with {h2 = }, {Is2[0] = }')
    print(f'Upper bound error for 1/3 rule: {E = }')
    print(f'Simpson\'s 1/3 error E = {abs(I0 - Is1[0])}')
    print(f'Simpson\'s 3/8 error E = {abs(I0 - Is2[0])}')

