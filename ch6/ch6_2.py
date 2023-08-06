import os
import sys
import matplotlib.pyplot as plt


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(DIR)


from la.linalg import *


def f1_diverge(args):
    return (10 - args[0] ** 2) / args[1]


def f2_diverge(args):
    return 57 - 3 * args[0] * args[1] ** 2


def f1(args):
    return (10 - args[0] * args[1]) ** 0.5


def f2(args):
    return ((57 - args[1]) / (3 * args[0])) ** 0.5



if __name__ == "__main__":

    # Example of nonlinear system of equations
    '''
        x1 ** 2 + x1 * x2 = 10
        x2 + 3 * x1 * x2 ** 2 = 57
    '''

    x0 = [1.5, 3.5]
    # Divergence functions
    # x = fixed_point([f1_diverge, f2_diverge], x0, 0.0001)

    # Convergence
    x = fixed_point([f1, f2], x0, 0.0001)
    print(str([round(xi, 3) for xi in x])[1:-1])
