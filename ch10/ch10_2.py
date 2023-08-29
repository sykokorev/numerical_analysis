import os
import sys
import math
import numpy as np


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(DIR)


from diff.naive_diff import *


def func(x, a0, a1, a2):
    return a0 * math.e ** (a1 * x) + a2


def func1(x):
    return 1.2 * x ** 3


if __name__ == "__main__":

    args = [0.2, 0.1, 0.3]
    x = 1.2

    d = diff(func, x, args=args)
    print(d)

    d2 = ndiff(func, x, 2, args)
    print(d2)

    d3 = ndiff(func, x, 3, args)
    print(d3)

    d = diff(func1, x)
    print(d)
    d2 = ndiff(func1, x, 2)
    print(d2)
    d3 = ndiff(func1, x, 3)
    print(d3)
