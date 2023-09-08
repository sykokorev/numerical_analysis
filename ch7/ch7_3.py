import os
import sys
import numpy as np
import matplotlib.pyplot as plt


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(DIR)


from interpolate.interpolate import *


if __name__ == "__main__":

    line = np.array([[1.0, 1.0], [2.0, 5.0]])
    second_order = np.array([[1.0, 1.0], [2.0, 5.0], [3.0, 2.0]])
    third_order = np.array([[1.0, 1.0], [2.0, 5.0], [3.0, 2.0], [3.2, 7.0]])
    forth_order = np.array([[1.0, 1.0], [2.0, 5.0], [3.0, 2.0], [3.2, 7.0], [3.9, 4.0]])
    xp = np.linspace(0, 4, 50)

    fig, ax = plt.subplots()
    ax.scatter(line[:,0], line[:,1], marker='o', color='r', label='Line')

    yp = []
    for xi in xp:
        yp.append(Lagrangeinterp(line[:,0], line[:,1], xi))
    ax.plot(xp, yp, marker='', color='r', label='Line')

    yp = []
    for xi in xp:
        yp.append(Lagrangeinterp(second_order[:,0], second_order[:,1], xi))
    ax.scatter(second_order[:,0], second_order[:,1], marker='o', color='g', label='2d order')
    ax.plot(xp, yp, marker='', color='g', label='2d order')

    yp = []
    for xi in xp:
        yp.append(Lagrangeinterp(third_order[:,0], third_order[:,1], xi))
    ax.scatter(third_order[:,0], third_order[:,1], marker='o', color='k', label='3d order')
    ax.plot(xp, yp, marker='', color='k', label='3d order')

    yp = []
    xp = np.linspace(0.8, 4, 50)
    for xi in xp:
        yp.append(Lagrangeinterp(forth_order[:,0], forth_order[:,1], xi))
    ax.scatter(forth_order[:,0], forth_order[:,1], marker='o', color='y', label='4th order')
    ax.plot(xp, yp, marker='', color='y', label='4th order')

    ax.grid(True)
    fig.legend()
    plt.show()
