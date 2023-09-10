import numpy as np

from integration.integrate import *
from la.linalg import *




if __name__ == "__main__":

    num_points = 3
    args = legendre_poly_coef(num_points)
    intervals = np.linspace(-1, 1, num_points + 1)
    print(intervals)
    xs = []
    for i, it in enumerate(intervals):
        try:
            xi = newton_raphson2(legendre_poly, it, args)
            xs.append(xi)
        except RuntimeError:
            print('Root no found')

    for xi in xs:
        print(xi)
