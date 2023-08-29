import numpy as np
import matplotlib.pyplot as plt
import math


from scipy.optimize import curve_fit


from interpolate.optimize import *


def cubic(x, a0, a1, a2, a3):
    return a0 + a1 * x + a2 * x ** 2 + a3 * x ** 3


def exp(x, a, b):
    return a * math.e ** (b * x)


def func(x, a, b, c):
    return a * np.exp(-b * x) + c


def mm_kinetic(x, b1, b2):
    return b1 * x / (b2 + x)


def Rosenbrock(xy, a, b):
    return (a - xy[0]) ** 2 + b * (xy[1] - xy[0] ** 2) ** 2


DIR = os.path.abspath('tmp')


if __name__ == "__main__":

    f = func
    params = [0.1, 0.5, 0.3]

    for i in range(10):

        print(f'Current iter {i}')
        
        xdata = np.linspace(0, 4, 50)
        y = f(xdata, *params)
        rng = np.random.default_rng()
        y_noise = 0.2 * rng.normal(size=xdata.size)
        ydata = y + y_noise
        
        try:
            popt, pcov = curve_fit(f, xdata, ydata)
            ysc = [f(xi, *popt) for xi in xdata]
        except RuntimeError:
            print('Runtime Error scipy', popt)
            popt = None

        try:
            popt_lm = lm(f, xdata, ydata)
            ylm = [f(xi, *popt_lm) for xi in xdata]
        except RuntimeError:
            print('Runtime Error lm')
            popt_lm = None

        if hasattr(popt, '__iter__') and hasattr(popt_lm, '__iter__'):
            ch = [(p - pn) > 1e-1 for p, pn in zip(popt, popt_lm)]

            df = sum([(f(xi, *popt_lm) - yi) ** 2 for xi, yi in zip(xdata, ydata)])
            df_sc = sum([(f(xi, *popt) - yi) ** 2 for xi, yi in zip(xdata, ydata)])

            print(df, df_sc, df - df_sc)
            fig, ax = plt.subplots()
            ax.scatter(xdata, ydata, marker='o', color='g', label='Noise data')
            ax.plot(xdata, ysc, marker='', color='r', label='scipy regression')
            ax.plot(xdata, ylm, marker='', color='k', label='my optimization')

            ax.grid(True)
            fig.legend()
            fig.savefig(os.path.join(DIR, f'{i}.png'))
            plt.close()
