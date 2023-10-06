import numpy as np


from de.ode import solve_bvp


def tempeq(x, y):
    t1 = y[1]
    t0 = y[0]
    return np.array([t1, 0.1 * (t0 - 20)])


def bcres(ya, yb):
    return np.array([ya[0] - 40, yb[0] - 200])


if __name__ == "__main__":

    x = np.linspace(0, 10, 6)
    y = np.zeros((2, x.size))
    print(y)

    sol = solve_bvp(tempeq, bcres, x, y)
    print(sol)
