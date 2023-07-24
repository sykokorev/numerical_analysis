import os
import sys


LIBDIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__),
    os.pardir)
)
sys.path.append(LIBDIR)

from la.linear_algebra import *



if __name__ == "__main__":

    pass
