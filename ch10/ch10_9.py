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


if __name__ == "__main__":

    args = [1.5, 2.5, 3.0, 0.15]

    x = np.linspace(0, 486, 20)
    y = [0.0 for i in range(x.size)]
    rng = np.random.default_rng()
    xdata = x + 5 * rng.normal(size=x.size)
    ydata = poly1D(xdata, *args)

    dxdy1 = deriv1D(xdata, ydata)

    fig, ax = plt.subplots()
    ax.plot(xdata, ydata, marker='o', color='k')

    plt.show()
    