import os
import sys
import matplotlib.pyplot as plt


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(DIR)


from la.linalg import *


def f1(args):
    return args[0] ** 2 + args[0] * args[1] - 10


def f2(args):
    return args[1] + 3 * args[0] * args[1] ** 2 - 57


def jacobian(args):
    return [
        [2 * args[0] + args[1], args[0]],
        [3 * args[1] ** 2, 1 + 6 * args[0] * args[1]]
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
    x = newton_raphson([f1, f2], jacobian, x0, 0.0001)
    print(x)
