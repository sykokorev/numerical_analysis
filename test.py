from __future__ import annotations

import os
import sys
import sympy

from sympy.core import mul, Expr


LIBDIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__),
    os.pardir)
)
sys.path.append(LIBDIR)

from symbols.symbol import *


if __name__ == "__main__":

    x = Symbol('x')
    z = Symbol('z')

    exp = x ** 2
    print(exp)
