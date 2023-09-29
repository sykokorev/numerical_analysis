import os
import sys
import matplotlib.pyplot as plt


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(DIR)


from la.linalg import *


def f1(x0, x1, x2):
    return x0 ** 2 + x1 * x0 - 10 * x2


def f2(x0, x1, x2):
    return x1 + 3 * x0 * x1 ** 2 - 57 * x2


def f3(x0, x1, x2):
    return x0 ** 3 + x1 * x0 - 2 * x2 - 12


def f4(x1, x2):
    return x1 + x2 - 1


def f5(x1, x2):
    return x1 ** 2 - x2 - 5


def jac(x0, x1):
    return [
        [2 * x0 + x1, x0],
        [3 * x0 ** 2 + x1, 1 + 6 * x0 * x1 ** 2]
    ]


if __name__ == "__main__":

    '''
    Solving nonlinear system of equations by using Newton-Raphson method
        x1 ** 2 + x1 * x2 = 10
        x2 + 3 * x1 * x2 ** 2 = 57
        {f}= {
            x1 ** 2 + x1 * x2 - 10
            x2 + 3 * x1 * x2 ** 2 - 57
        }
        Jacobian [J] = [
            2 * x1 + x2,    x1
            3 * x2 ** 2,    1 + 6 * x1 * x2
        ]
    '''

    x0 = [1.5, 3.5]
    x = newton_raphson([f1, f2, f3])
    print(f1(*x), f2(*x), f3(*x))
    print(x)
    x = newton_raphson([f1])
    print(x)
