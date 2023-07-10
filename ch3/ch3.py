import numpy as np
import matplotlib.pyplot as plt
import random
import os
import sys


LIBDIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__),
    os.pardir)
)
sys.path.append(LIBDIR)


from lib.consts import COLORS, MARKERS
from lib.classes import *
from lib.functions import *


if __name__ == "__main__":

    # Exponantial function. Errors analysis.
    func = Exponential(0.5, 1)
    n = np.linspace(1, 50, 50)
    
    errors = []
    y = []
    x0 = 8
    x1 = 5

    for i in n:
        y5 = taylor_series(func, x1, x0, i)
        y.append(y5)
        errors.append(max_error(func, x1, x0, i))

    fig, axs = plt.subplots(2, 1)
    axs[0].plot(n, y, linestyle='-', marker='*', color='red')
    axs[0].set_title('Function approxiamtion')
    axs[1].plot(n, errors, linestyle='-', marker='o', color='blue')
    axs[1].set_title('Residuals')
    axs[0].grid(True)
    axs[1].grid(True)

    plt.show()


    # Polynomial function. Approximation analysis
    poly = Polynomial([1.2, -0.25, -0.5, -0.15, -0.1])

    fig, ax = plt.subplots()

    x = np.linspace(-2, 2, 21)
    x0 = 0
    n = np.linspace(0, 4, 5)
    y = []

    for xi in x:
        y.append(poly(xi))
    ax.plot(x, y, label='exact', color='k', marker='o')

    for i in n:
        y = []
        for xi in x:
            y.append(taylor_series(poly, xi, x0, i))
        color = random.randint(0, len(COLORS)-1)
        marker = random.randint(0, len(MARKERS) - 1)
        ax.plot(
            x, y, color=COLORS[color], 
            marker=MARKERS[marker], label=f'{int(i)}th order'
        )
        ax

    ax.grid()
    ax.set_xlabel('x')
    ax.set_ylabel('Approximation to f(x) = -0.1$x^4$ - 0.15$x^3$ - 0.5$x^2$ - 0.25x + 1.2')
    fig.legend()
    plt.show()    
