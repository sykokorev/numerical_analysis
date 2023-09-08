from __future__ import annotations

import os
import sys
import numpy as np

from sympy.core import mul, Expr


LIBDIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__),
    os.pardir)
)
sys.path.append(LIBDIR)

from integration.integrate import romberg


if __name__ == "__main__":

    gaussian = lambda x: 1 / np.sqrt(np.pi) * np.exp(-x ** 2)

    I = romberg(gaussian, 0, 1)
    print(I)

    
