import os
import sys

import numpy as np
import scipy as sp


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
    # print([round(xi, 3) for xi in x])
    # print(np.linalg.solve(a, b))

    a = [
        [2, 1, 1],
        [5, 2, 2],
        [4, 3, 2]
    ]

    b = [5, 6, 3]

    print('Initial matrix')
    for ai in a:
        print(ai)

    lower, upper = lu(a)

    print('Upper triangular matrix')
    for ui in upper:
        print(ui)

    print('Lower triangular matrix')
    for li in lower:
        print(li)
    
    print('Checking LU by multiply L * U')
    init = mat_mul(lower, upper)
    for i in init:
        print(i)

    # Solving equation
    x = solveLU(a, b)
    for i, xi in enumerate(x, 1):
        print(f'x{i} = {xi}', end='\t')
    print()

    # Cholesky Factorization
    a = [
        [1, 3, 2],
        [3, 13, 8],
        [2, 8, 6]
    ]
    b = [2, 4, 1]

    x = cholab(a, b)
    for i, xi in enumerate(x, 1):
        print(f'x{i} = {xi}', end='\t')
    print()
