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

    data = np.array([
        [1.0, 1.93], [1.1, 1.61], [1.2, 2.27],
        [1.3, 3.19], [1.4, 3.19], [1.5, 3.71],
        [1.6, 4.29], [1.7, 4.95], [1.8, 6.07],
        [1.9, 7.48], [2.0, 8.72], [2.1, 9.34],
        [2.2, 11.62]
    ])
    x, y = data[:, 0], data[:, 1]

    # Exponential model
    (a, b), r = exp(x, y)

    print('The parameters of model y(x) = b * e ** (a * x) for data:')
    for xi, yi in zip(x, y):
        print(f'\t{(xi, yi)}')
    print(f'{round(a, 4) = }, {round(b, 4) = }, y(x) = {round(a, 4)} * e ** ({round(b, 4)} * x)')
    print(f'with the correlation coefficient {round(r, 3)}')

    xm = np.linspace(1.0, 2.2, 100)
    ym = []

    for xi in xm:
        ym.append(exp_model(a, b, xi))

    # Power model
    (a, b), r = power(x, y)

    ypm = []
    for xi in xm:
        ypm.append(power_model(a, b, xi))

    print('The parameters of model y(x) = b * x ** a for data:')
    for xi, yi in zip(x, y):
        print(f'\t{(xi, yi)}')
    print(f'{round(a, 4) = }, {round(b, 4) = }, y(x) = {round(b, 4)} * x ** {round(a, 4)})')
    print(f'with the correlation coefficient {round(r, 3)}')

    fig, ax = plt.subplots()
    ax.scatter(x, y, marker='o', color='g', label='Observed data')
    ax.plot(xm, ym, marker='', color='k', label='Modeled data $y(x) = b \star e^{a \star x}$')
    ax.plot(xm, ypm, marker='', color='b', label='Modeled data $y(x) = b \star x^{a}$')

    ax.grid(True)
    fig.legend()

    plt.show()
