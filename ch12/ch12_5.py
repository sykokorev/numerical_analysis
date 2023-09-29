import os
import sys
import numpy as np
import matplotlib.pyplot as plt


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__), os.pardir
))
sys.path.append(DIR)


from de.ode import solve_ipv


def hc(x):
    return -20.7302 * np.sinh(0.316228 * x) + 20 * np.cosh(0.316288 * x) + 20


def hc2(x):
    return -4.80593 * np.sinh(0.316228 * x) + 20 * np.cosh(0.316288 * x) + 20


def y2p(x, y):
    return [y[1], 2 * x + x * y[0] - x * y[1]]


def heat_conserv(x, t, h, ta):
    '''
        The law of heat conservation T''(x) + h (Ta - T) = 0
            x is a length of 1D rod (the independent variable)
            t is a temperature of 1D rod (the dependent variable)
        t : list
            t[0] = T(x)
            t[1] = T'(x)
    '''
    return [t[1], h * (t[0] - ta)]


if __name__ == "__main__":

    h = 0.1
    ta = 20
    ic = [40, -1.5]
    xspan = np.linspace(0, 10, 5)
    xsol, tsol = solve_ipv(heat_conserv, [xspan[0], xspan[-1]], ic, teval=xspan, args=[h, ta])

    xs = np.linspace(0, 10, 100)
    ts = hc2(xs)
    print(tsol)


    fig, ax = plt.subplots()
    ax.scatter(xsol, tsol[0], marker='o', color='b', label='solve_ivp')
    ax.plot(xs, ts, 'r--', label='Exact solution')

    ax.grid(True)

    plt.legend()
    plt.show()
