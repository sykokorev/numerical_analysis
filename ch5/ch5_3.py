import os
import sys


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(DIR)


from la.linalg import *


if __name__ == "__main__":

    a = [
        [3.0, -0.1, -0.2],
        [0.1, 7.0, -0.3],
        [0.3, -0.2, 10]
    ]
    b = [7.85, -19.3, 71.4]

    # Jacobi method
    x, error = jacobi(a, b, [1.0, 1.0, 1.0], 0.0001)

    for i, xi in enumerate(x, 1):
        print(f'x{i} = {round(xi, 4)}', end='\t')
    print(f'Relative error {error}')

    # Gauss-Seidel method
    x, error = gauss_seidel(a, b, [1.0, 1.0, 1.0], 0.0001)

    for i, xi in enumerate(x, 1):
        print(f'x{i} = {round(xi, 4)}', end='\t')
    print(f'Relative error {error}')
