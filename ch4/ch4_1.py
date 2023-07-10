from __future__ import annotations

import os
import sys


LIBDIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__),
    os.pardir)
)
sys.path.append(LIBDIR)

from lib.classes import *
from lib.functions import bisection


if __name__ == "__main__":

    a = [-0.6, -0.3, 0.6]
    b = [-0.5, -0.2, 0.7]
    e = 0.0005
    cos = Cos(2)
    sin = Sin(5)

    func = Addition([sin, cos])

    print(f'The roots of the equation {func} with maximal error {e} are:')

    for ai, bi in zip(a, b):
        print(f'\tInside the interval [{ai, bi}]: {bisection(func, ai, bi, e=e)}')

    a = [1]
    b = [2]
    func = Polynomial([-10, 0, 2, 1])
    print(f'The roots of the equation: {func} are:')
    for ai, bi in zip(a, b):
        print(f'\tInside the interval [{ai, bi}]: {bisection(func, ai, bi, e=e)}')
