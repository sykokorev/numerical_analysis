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


from integration.integrate import gauss_quad



def func1(x):
    return 2000 * math.log(140000 / (140000 - 2100 * x)) - 9.8 * x


def func2(x):
    return 1 / (2 + x)


def func2int(x):
    return math.log(x)


if __name__ == "__main__":

    I = func2int(3) - func2int(1)
    print(f'Exact integral I = {I}')

    for i in range(2, 40):

        I1 = gauss_quad(f=func1, a=8, b=30, n=i)
        I2 = gauss_quad(f=func2, n=i)
        print(i, I1)
        print(i, I2, abs(I2 - I) / I)
