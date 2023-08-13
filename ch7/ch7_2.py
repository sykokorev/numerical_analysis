import os
import sys
import numpy as np
import matplotlib.pyplot as plt


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(DIR)


from interpolate.interpolation import *


if __name__ == "__main__":

    points = np.array([
        [1.0, 1.0],
        [2.0, 5.0],
        [3.0, 2.0], 
        [3.2, 7.0],
        [3.9, 4.0]
    ])
    xp = np.linspace(0.8, 4, 100)
    yp = []

    for xi in xp:
        yp.append(Newtoninterp(points[:, 0], points[:, 1], xi))

    fig, ax = plt.subplots()
    ax.scatter(points[:, 0], points[:, 1], marker='o', color='g')
    ax.plot(xp, yp, marker='', color='k')

    ax.grid(True)
    plt.show()
