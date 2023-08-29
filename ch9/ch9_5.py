import os
import sys
import numpy as np
import matplotlib.pyplot as plt


from scipy.optimize import curve_fit


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(DIR)


from interpolate.optimize import *


def fit_f(x, a, b, c, d):
    return a*(x-b)**2+c+d*0.0001*np.cos(x)


def cubic(x, a0, a1, a2, a3):
    return a0 + a1 * x + a2 * x ** 2 + a3 * x **3


if __name__ == "__main__":

    xdata = np.array([ 0.23547456, 0.15789474, 0.31578947, 0.47368421, 0.63157895, 
                    0.78947368, 0.94736842, 1.10526316, 1.26315789, 1.42105263, 
                    1.57894737, 1.73684211, 1.89473684, 2.05263158, 2.21052632, 
                    2.36842105, 2.52631579, 2.68421053, 2.84210526, 3.45454545 ])
    ydata = np.array([ 2.95258285, 2.49719803, -2.1984975, -4.88744346, -7.41326345, 
                    -8.44574157, -10.01878504, -13.83743553, -12.91548145, -15.41149046, 
                    -14.93516299, -13.42514157, -14.12110495, -17.6412464 , -16.1275509 , 
                    -16.11533771, -15.66076021, -13.48938865, -11.33918701, -11.70467566])

    args = [1.0, 1.0, 1.0, 1.0]

    popt = lm(fit_f, xdata, ydata)
    popt_sc, pcov = curve_fit(fit_f, xdata, ydata)

    y = []
    ysc = []

    for x in xdata:
        y.append(fit_f(x, *popt))
        ysc.append(fit_f(x, *popt_sc))
    

    fig, ax = plt.subplots()
    ax.scatter(xdata, ydata, marker='o', color='g', label='Observed data')
    ax.plot(xdata, y, marker='', color='r', label='Approximated data')
    ax.plot(xdata, ysc, marker='', color='b', label='Approximated data scipy')

    ax.grid(True)
    fig.legend()

    plt.show()
