import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt


from scipy.integrate import romberg as romb


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(DIR)

from integration.integrate import romberg, trapezoid


def exp(x):
    return np.exp(x)


def exp2(x):
    return 5 * x * np.exp(-2 * x)


def func(x):
    return 2 + 2 * x + x ** 2 + math.sin(2 * math.pi * x) + math.cos(4 * math.pi * x)


def funcint(x):
    return 2 * x + x ** 2 + (x ** 3) / 3 - math.cos(2 * math.pi * x) / (2 * math.pi) + math.sin(4 * math.pi * x) / (4 * math.pi)


if __name__ == "__main__":

    a = 0.0
    b = 1.5

    Iex = funcint(b) - funcint(a)
    print(f'Exact itegration: {round(Iex, 5)}')

    n = [1, 2, 4, 8, 16, 32, 64, 128]
    for ni in n:
        I = trapezoid(func, a, b, ni)
        es = abs(I - Iex)
        er = es / Iex
        print(ni, round(I, 5), round(er, 7))

    Ir = romberg(func, a, b, rtol=0.00042)
    print(f'Romberg rule: {Ir}')
    print(f'Romberg rule with scipy lib: {romb(func, a, b, rtol=0.00042)}')
