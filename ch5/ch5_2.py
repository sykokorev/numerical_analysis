import os
import sys

import numpy as np


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(DIR)

from la.linalg import *


if __name__ == "__main__":

    # Cramer's Rule
    a = [
        [1, 2],
        [-2, 1]
    ]
    b = [-1, 2]

    x = cramer(a, b)
    # print(x)
    a = [
        [1, 2, 3],
        [2, 1, 4],
        [1, 3, 5]
    ]
    b = [-4, 8, 0]
    x = cramer(a, b)
    # print(x)


    # Naive-Gaussian
    a = [
        [4, 1, 2],
        [2, 3, 4],
        [-1, 2, 6]
    ]
    b = [6, 4, 8]

    x = naive_gaussian(a, b)

    # Gauss-Jordan Elimination
    x = gauss_jordan(a, b)
    print([round(xi, 3) for xi in x])
    print(np.linalg.solve(a, b))