import os
import sys
import math
import numpy as np


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(DIR)


from diff.diff import *


def func(x, a0, a1, a2):
    return a0 * math.e ** (a1 * x) + a2


def func1(x):
    return 1.2 * x ** 3


if __name__ == "__main__":

    x = 1.2
    args = [0.2, 0.1, 0.3]

    d1 = diff(func, x, 1, args)
    d2 = diff(func, x, 2, args)
    d3 = diff(func, x, 3, args)
    d4 = diff(func, x, 4, args)
    print(d1)
    print(d2)
    print(d3)
    print(d4)
