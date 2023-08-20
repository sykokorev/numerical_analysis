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

    y0 = runge(1, 1, 5, 0.0)
    y1 = runge(1, 1, 5, -0.01)
    y2 = runge(1, 1, 5, -0.001)
    y3 = runge(1, 1, 5, 0.001)
    y4 = runge(1, 1, 5, 0.01)

    points = np.array([
        [-1.0, 0.038], [-0.8, 0.058],
        [-0.60, 0.1], [-0.4, 0.2],
        [-0.2, 0.5], [0.0, y0],
        [0.2, 0.5], [0.4, 0.2],
        [0.6, 0.1], [0.8, 0.058],
        [1.0, 0.038]
    ])

    # points = np.array([
    #     [-1.0, 0.038], [-0.8, 0.058],
    #     [-0.6, 0.1], [-0.4, 0.2],
    #     [-0.2, 0.5]
    # ])

    # interp = Interpolate(points, 1)
    # interp.eval()

    interp = Interpolate(points, 2)
    interp.eval()

    print(interp.intervals)

    xp = np.linspace(-1.0, 1.0, 1000)
    yp = []
    y1 = []
    for xi in xp:
        yp.append(runge(1, 1, 5, xi))
        y1.append(interp(xi)[1])

    xp = np.linspace(-1.0, 1.0, 1000)
    yp = []
    y1 = []
    for xi in xp:
        yp.append(runge(1, 1, 5, xi))
        y1.append(interp(xi)[1])
    
    print('Spline consists from:')
    print(interp)

    fig, ax = plt.subplots()
    ax.scatter(points[:, 0], points[:, 1], marker='o', color='g', label='Points')
    ax.plot(xp, y1, marker='', color='k', label=f'Polynomial spline degree {interp.n}')

    ax.grid(True)
    fig.legend()
    plt.show()
