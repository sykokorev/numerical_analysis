import os
import sys
import matplotlib.pyplot as plt
import numpy as np

DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(DIR)


from interpolate.interpolation import *


def polynomial(coef: list, x: list):
    return [sum([coef[j] * xi ** j for j in range(len(coef))]) for xi in x]



if __name__ == "__main__":

    points = [(1.1, 1.0), (2.0, 2.1), (3.0, 5.0), (3.4, 7.0)]

    xp = [xi[0] for xi in points]
    yp = [yi[1] for yi in points]

    a = polyfit(points)
    print(a)
    x = np.linspace(0, 5, 1000)
    y = polynomial(a, x)

    fig, ax = plt.subplots()
    ax.scatter(xp, yp, color='g', marker='o')
    ax.plot(x, y, linestyle='-', marker='', color='k')

    ax.grid(True)
    plt.show()
