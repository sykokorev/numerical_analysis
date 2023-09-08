import os
import sys


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))

sys.path.append(DIR)


from la.linear_algebra import zeros


def richardson(f: callable, arg: float | int, h: float, args: list = [], k: int = 2):
    
    '''
        Richardson's extrapolation
        --------------------------
            f: callable: function need to be extrapolated
            arg: int | float: function argument
            h: float: initial step
            args: list (default []): functio arguments (default there are no arguments. [])
            k: int (default 2): order of error terms (O(h ** k))
        --------------------------
        return tuple. Returns tuple of list of functions value, relative error
    '''

    r = zeros(k, axis=2)
    for j in range(k):
        for i in range(k):

            if j == 0:
                x = arg + h if not i else arg + h / (2 ** i)
                r[i][j] = f(x, *args)
            else:
                const = 4 ** j
                r[i][j] = (const * r[i + 1][j - 1] - r[i][j - 1]) / (const - 1)
    
    tol = abs(r[0][k - 1] - r[0][k]) / r[0][k - 1]

    return r, tol
