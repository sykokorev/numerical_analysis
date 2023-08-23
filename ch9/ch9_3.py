import os
import sys
import numpy as np
import matplotlib.pyplot as plt


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(DIR)


from interpolate.curve_fit import *


def f1(*args):
    return 1


def f2(x):
    return x


def f3(x):
    return np.cos(np.pi * x)


def c1(*args):
    return 1


def c2(x):
    return x


def c3(x):
    return x ** 2


def c4(x):
    return x ** 3


if __name__ == "__main__":

    data = np.array([
        [1.0, 0.5], [2.0, 2.5], [3.0, 2.0],
        [4.0, 4.0], [5.0, 3.5],
        [6.0, 6.0], [7.0, 5.5]
    ])

    func = [f1, f2, f3]
    x = data[:, 0]
    y = data[:, 1]

    (a1, a2, a3), r = ext_linear(func, x, y)

    print(f'Parameters of the model y(x) = a1 + a2 * x + a3 * cos(pi * x) for data set:')
    for xi, yi in zip(x, y):
        print(f'\t{(xi, yi)}')
    print(f'a1 = {round(a1, 5)}, a2 = {round(a2, 5)}, a3 = {round(a3, 5)}')
    print(f'Correlation coefficient R squared = {round(r, 5)}', end='\n\n')

    xm = np.linspace(1.0, 7.0, 10000)
    ym = []
    for xi in xm:
        ym.append(linear_model(xi, [f1, f2, f3], [a1, a2, a3]))

    fig, ax = plt.subplots()
    ax.scatter(x, y, marker='o', color='g', label='Observed data')
    ax.plot(xm, ym, marker='', color='k', label='Modeled data')

    ax.grid(True)
    fig.legend()

    plt.show()


    # Cubic polynomial model
    data = np.array([
        [1.0, 1.93], [1.1, 1.61], [1.2, 2.27], [1.3, 3.19],
        [1.4, 3.19], [1.5, 3.71], [1.6, 4.29], [1.7, 4.95], 
        [1.8, 6.07], [1.9, 7.48], [2.0, 8.72], [2.1, 9.34], 
        [2.2, 11.62]
    ])

    x = data[:, 0]
    y = data[:, 1]
    f = [c1, c2, c3, c4]

    (a1, a2, a3, a4), r = ext_linear(f, x, y)

    print(f'Parameters of the model y(x) = a1 + a2 * x + a3 * x ** 2 + a4 * x ** 3 for data set:')
    for xi, yi in zip(x, y):
        print(f'\t{(xi, yi)}')
    print(f'a1 = {round(a1, 5)}, a2 = {round(a2, 5)}, a3 = {round(a3, 5)}, a4 = {round(a4, 5)}')
    print(f'Correlation coefficient R squared = {round(r, 5)}')

    xm = np.linspace(1.0, 2.2)
    ym = []
    for xi in xm:
        ym.append(linear_model(xi, f, [a1, a2, a3, a4]))


    fig, ax = plt.subplots()
    ax.scatter(x, y, marker='o', color='g', label='Observed data')
    ax.plot(xm, ym, marker='', color='k', label='Modeled data')

    ax.grid(True)
    fig.legend()

    plt.show()    
