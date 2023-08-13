import os
import sys
import numpy as np
import matplotlib.pyplot as plt


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(DIR)


from interpolate.piecewice import *


def runge(a, b, c, x):
    return (a / (b + (c * x) ** 2))


if __name__ == "__main__":

    # points = np.array([
    #     [-1.0, 0.038], [-0.8, 0.058],
    #     [-0.60, 0.1], [-0.4, 0.2],
    #     [-0.2, 0.5], [0.0, 1.0],
    #     [0.2, 0.5], [0.4, 0.2],
    #     [0.6, 0.1], [0.8, 0.058],
    #     [1.0, 0.038]
    # ])

    points = np.array([[-1.0, 0.038], [-0.8, 0.058], [-0.6, 0.1], [-0.4, 0.2]])

    interp = Interpolate(points, 1)
    interp.eval()

    xp = np.linspace(-1.0, 1.0, 1000)
    yp = []
    y1 = []
    for xi in xp:
        yp.append(runge(1, 1, 5, xi))
        y1.append(interp(xi)[1])

    interp.n = 2

    # fig, ax = plt.subplots()

    # ax.plot(xp, yp, marker='', color='k', label='Runge function')
    # ax.plot(xp, y1, marker='', color='b', label='Piecewice Linear')
    # ax.grid(True)
    # fig.legend()

    # plt.xticks(np.arange(min(xp), max(xp)+0.2, 0.2))
    # plt.yticks(np.arange(min(yp), max(yp)+0.1, 0.1))
    # plt.show()

