import numpy as np
import matplotlib.pyplot as plt


from scipy.optimize import curve_fit


from la.optimize import *


def cubic(x, a0, a1, a2, a3):
    return a0 + a1 * x + a2 * x ** 2 + a3 * x ** 3


def func(x, a, b, c):
    return a * np.exp(-b * x) + c


def jac(x, a0, a1, a2, a3):
    return [1, x, x ** 2, x ** 3]


if __name__ == "__main__":

    # for i in range(1000):
        # print(i)
        xdata = np.linspace(0, 4, 50)
        y = func(xdata, 15, 25, 31)
        rng = np.random.default_rng()
        y_noise = 0.2 * rng.normal(size=xdata.size)
        ydata = y + y_noise

        popt, pcov = curve_fit(func, xdata, ydata)
        ysc = [func(xi, *popt) for xi in xdata]
        print(f'{popt=}')

        popt_lm, error = lm(func, xdata, ydata)
        ylm = [func(xi, *popt_lm) for xi in xdata]
        print(f'{popt_lm=}')
        print([p - pn for p, pn in zip(popt, popt_lm)])

        fig, ax = plt.subplots()
        ax.scatter(xdata, ydata, marker='o', color='g', label='Noise data')
        ax.plot(xdata, ysc, marker='', color='r', label='scipy regression')
        ax.plot(xdata, ylm, marker='', color='k', label='my optimization')

        ax.grid(True)
        fig.legend()

        plt.show()

    # Test cubic function
    # xdata = np.linspace(0, 4, 50)
    # y = cubic(xdata, 50, 25, 51, 0.01)
    # rng = np.random.default_rng()
    # y_noise = rng.normal(scale=100, size=xdata.size)
    # ydata = y + y_noise

    # popt, error = lm(cubic, xdata, ydata)
    # print('Error: ', error)
    # popt_sc, pcov = curve_fit(cubic, xdata, ydata)
    # print(f'{popt_sc=}')
    # print(f'Optimal parameters of model function are a0 = {popt[0]}, a1 = {popt[1]}, a2 = {popt[2]}, a3 = {popt[3]}')
    # print(f'Model function y(x) = a0 + a1 * x + a2 * x3 ** 2 + a3 * x3 ** 3')

    # ylm = [cubic(xi, *popt) for xi in xdata]
    # ysc = [cubic(xi, *popt_sc) for xi in xdata]

    # fig, ax = plt.subplots()
    # ax.scatter(xdata, ydata, marker='o', color='g', label='Noise data')
    # ax.plot(xdata, ysc, marker='', color='r', label='scipy regression')
    # ax.plot(xdata, ylm, marker='', color='k', label='my optimization')

    # ax.grid(True)
    # fig.legend()

    # plt.show()

