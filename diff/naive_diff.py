def fact(n):

    if n == 0:
        return 1
    elif n < 0:
        return -1
    else:
        return n * fact(n - 1)


def binom(i: int, n: int):
    return fact(n) / (fact(i) * fact(n - i))


def fdiff(f, x: float, args: list = [], dx: float = 1e-3):
    '''
        Compute forward finite difference
            f: callable: The function to be differentiated
            x: float: point at which function to be differentiated
            args: list (default []): List of function's arguments
            dx: float (default 1e-8): delta x
    '''
    dx *= x
    return (f(x + dx, *args) - f(x, *args)) / dx


def bdiff(f, x: float, args: list = [], dx: float = 1e-3):
    '''
        Compute forward finite difference
            f: callable: The function to be differentiated
            x: float: point at which function to be differentiated
            args: list (default []): List of function's arguments
            dx: float (default 1e-8): delta x
    '''
    dx *= x
    return (f(x, *args) - f(x - dx, *args)) / dx


def diff(f, x: float, args: list = [], dx: float = 1e-3):
    '''
        Compute centered finite difference
            f: callable: The function to be differentiated
            x: float: point at which function to be differentiated
            args: list (default []): List of function's arguments
            dx: float (default 1e-8): delta x
    '''
    dx *= x
    return (f(x + dx, *args) - f(x - dx, *args)) / (2 * dx)


def ndiff(f, x: float, n: int = 1, args: list = [], dx: float = 1e-3):
    '''
        Compute n order derivative
            f: callable: The function to be differentiated
            x: float: point at which function to be differentiated
            n: order of derivatives
            args: list (default []): List of function's arguments
            dx: float (default 1e-4): delta x
    '''
    dx *= x
    return (1 / dx ** n) * sum([(-1) ** (k + n) * binom(k, n) * f(x + k * dx, *args) for k in range(n + 1)])
