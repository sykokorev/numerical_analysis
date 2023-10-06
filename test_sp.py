import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import solve_bvp



def tempext(x):
    return -4.80593 * np.sinh(0.316228 * x) + 20 * np.cosh(0.316288 * x) + 20


def ext2(x):
    return -20.7302 * np.sinh(0.316228 * x) + 20 * np.cosh(0.316288 * x) + 20


def tempeq(x, y):
    t1 = y[1]
    t0 = y[0]
    print(x, y)
    return np.array([t1, 0.1 * (t0 - 20)])


def bcres(ya, yb):
    print(ya, yb)
    print(f'bcres: {[ya[0] - 40, yb[0] - 200] = }')
    return np.array([ya[0] - 40, yb[0] - 200])


def bcres1(ya, yb):
    print(f'bcres1: {[ya[0] - 40, yb[1] + 3] = }')
    return np.array([ya[0] - 40, yb[1] + 3])


if __name__ == "__main__":

    x = np.linspace(0, 10, 5)
    xext = np.linspace(0, 10, 500)
    y = np.zeros((2, x.size))

    res = solve_bvp(tempeq, bcres, x, y)
    ext = tempext(xext)

    print(t := res.sol(x)[0])
    print(res.sol(x)[1])

    fig, ax = plt.subplots()
    ax.scatter(x, t)
    ax.plot(xext, ext)

    plt.show()

    x = np.linspace(0, 10, 19)
    y = np.zeros((2, x.size))
    t2 = solve_bvp(tempeq, bcres1, x, y)

    print(t := t2.sol(x)[0])

    fig, ax = plt.subplots()
    ax.plot(xext, ext2(xext))
    ax.scatter(x, t)

    plt.show()
