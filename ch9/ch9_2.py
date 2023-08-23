import os
import sys
import numpy as np
import matplotlib.pyplot as plt


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(DIR)


from interpolate.curve_fit import *


if __name__ == "__main__":

    points = np.array([
        [1, 0.5], [2, 2.5], [3, 2.0],
        [4, 4.0], [5.0, 3.5],
        [6, 6.0], [7, 5.5]
    ])

    (a, b), r = simple_linear(points[:, 0], points[:, 1])

    print(f'Parameters of the linear model are a = {a}, b = {b}. Correlation coefficient R^2 = {r}')
    
    y = []
    for x in points[:, 0]:
        y.append(simple_linear_model(a, b, x))

    fig, ax = plt.subplots()
    ax.scatter(points[:, 0], points[:, 1], marker='o', color='r', label='observated data')
    ax.plot(points[:, 0], y, marker='', color='k', label='linear model')
    
    ax.grid(True)
    fig.legend()

    plt.show()
