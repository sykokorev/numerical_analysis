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


from diff.deriv import *


def poly1D(x, a0, a1, a2, a3):
    return a0 + a1 * x + a2 * x ** 2 + a3 * x ** 3


def dpoly1D(x, a0, a1, a2, a3):
    return a1 + 2 * a2 * x + 3 * a3 * x ** 2


if __name__ == "__main__":

    args = [0.5, 2.5, 1.0, 2.5]

    x = np.linspace(0, 2, 50)
    y = [0.0 for i in range(x.size)]
    rng = np.random.default_rng()
    xdata = x + rng.normal(size=x.size, scale=0.01)
    ydata = poly1D(xdata, *args)

    dxdy1 = deriv1D(xdata, ydata)
    dxdy2 = dpoly1D(x, *args)

    fig, ax = plt.subplots()
    ax.plot(xdata, ydata, marker='o', color='k')
    ax.plot(xdata, dxdy1, marker='*', color='r')
    ax.plot(x, dxdy2, marker='*', color='g')

    plt.show()
    