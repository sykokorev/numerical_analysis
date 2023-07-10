from __future__ import annotations

import os
import sys


LIBDIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__),
    os.pardir)
)
sys.path.append(LIBDIR)

from lib.classes import *


if __name__ == "__main__":

    cos = Cos(1.0, 5.0)
    sin = Sin(1.0, 2.0)
    cos1 = Cos(2.0, 3.0)

    sin *= cos
    sin *= cos1
    print(sin)
    # print()
