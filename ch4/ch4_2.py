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

    pass
