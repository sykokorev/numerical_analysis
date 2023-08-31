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


from interpolate.optimize import lm
from diff.diff import *


def runge(x, a, b, c):
    return a / (b + (x * c) ** 2)


def poly3(x, a0, a1, a2, a3):
    return a0 + a1 * x + a2 * x ** 2 + a3 * x ** 3


def poly4(x, a0, a1, a2, a3, a4):
    return a0 + a1 * x + a2 * x ** 2 + a3 * x ** 3 + a4 * x ** 4


def poly5(x, a0, a1, a2, a3, a4, a5):
    return a0 + a1 * x + a2 * x ** 2 + a3 * x ** 3 + a4 * x ** 4 + a5 * x ** 5


def poly6(x, a0, a1, a2, a3, a4, a5, a6):
    return a0 + a1 * x + a2 * x ** 2 + a3 * x ** 3 + a4 * x ** 4 + a5 * x ** 5 + a6 * x ** 6


def poly10(x, a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10):
    return a0 + a1 * x + a2 * x ** 2 + a3 * x ** 3 + a4 * x ** 4 + a5 * x ** 5 + a6 * x ** 6 + a7 * x ** 7+ a8 * x ** 8 + a9 * x ** 9 + a10 * x ** 10


if __name__ == "__main__":

    # Runge function
    args = [1, 1, 5]
    x = np.linspace(-1, 1, 10001)
    xdata = np.linspace(-1, 1, 51)

    y = runge(x, *args)
    yd = runge(xdata, *args)
    ydata = runge(xdata, *args)

    popt, pcov = curve_fit(poly5, xdata, ydata)
    popt_lm = lm(poly5, xdata, ydata)

    ysc = poly5(xdata, *popt)
    ylm = poly5(xdata, *popt_lm)

    d1 = []
    d2 = []
    d3 = []
    d4 = []
    for xi in x:
        d1.append(derivative(runge, xi, 1, args))
        d2.append(derivative(runge, xi, 2, args))
        d3.append(derivative(runge, xi, 3, args))
        d4.append(derivative(runge, xi, 4, args))

    x1, y1 = differnce(xdata, yd, 1)
    x2, y2 = differnce(xdata, yd, 2)
    x3, y3 = differnce(xdata, yd, 3)
    x4, y4 = differnce(xdata, yd, 4)

    fig, axs = plt.subplots(nrows=3, ncols=2)
    axs[0][0].plot(x, y, marker='', color='b', label='Runge func')
    axs[0][0].plot(xdata, ysc, marker='*', color='r', label='scipy')
    axs[0][0].plot(xdata, ylm, marker='>', color='k', label='LM')


    axs[0][0].scatter(xdata, ydata, marker='o', label='Runge func')
    axs[1][0].plot(x, d1, marker='', color='r', label='1st derivative')
    axs[1][0].scatter(x1, y1, marker='o', color='k', label='1st derivative app')
    axs[0][1].plot(x, d2, marker='', color='r', label='2d derivative')
    axs[0][1].scatter(x2, y2, marker='o', color='k', label='2d derivative app')
    axs[1][1].plot(x, d3, marker='', color='r', label='3d derivative')
    axs[1][1].scatter(x3, y3, marker='o', color='k', label='3d derivative app')
    axs[2][0].plot(x, d4, marker='', color='r', label='3d derivative')
    axs[2][0].scatter(x3, y4, marker='o', color='k', label='4d derivative app')

    axs[0][0].grid(True)
    axs[0][0].set_title('Runge function')
    axs[1][0].grid(True)
    axs[1][0].set_title('1st derivative')
    axs[0][1].grid(True)
    axs[0][1].set_title('2d derivative')
    axs[1][1].grid(True)
    axs[1][1].set_title('3d detivative')
    axs[2][0].grid(True)
    axs[2][0].set_title('4d detivative')

    fig.legend()
    plt.show()


