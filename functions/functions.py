import os
import sys
import numpy as np


def fact(n: int):
    
    if n <= 1:
        return 1
    return n * fact(n-1)


def cos(x: float, n: int) -> float:

    if n < 1:
        return 1

    return ((-1) ** n) * (x ** (2 * n)) / fact(2 * n) + cos(x, n - 1)


def sin(x: float, n: int) -> float:

    if n < 1:
        return x

    return ((-1) ** n) * (x ** (2 * n + 1)) / fact(2 * n + 1) + sin(x, n - 1)


def exp(x: float, n: int) -> float:

    if n < 1:
        return 1
    
    return (x ** n) / fact(n) + exp(x, n - 1)


def differential(f, x: float, n: int):

    if n <= 0:
        return f(x)

    return differential(f.derivative(), x, n - 1)


def taylor_series(f, x: float, x0: float, n: int):

    if n == 0:
        return differential(f, x0, n)
   
    return (differential(f, x0, n) * (x - x0) ** n) / fact(n) + taylor_series(f, x, x0, n - 1)


def max_error(f, x: float, x0: float, n: int) -> list:

    return (max(
        [abs(differential(f, x, n + 1)), 
         abs(differential(f, x0, n + 1))]) * (x - x0) ** (n + 1)) / fact(n + 1)


def bisection(f, a: float, b: float, iter: int = 1000, e: float = 10 ** -6):

    it = iter - 1
    if it == 0:
        raise RecursionError("The maximum numbers of iterations has been reached. No roots found.")

    if np.sign(f(a)) == np.sign(f(b)):
        raise Warning("There are no roots that can be found by the bisection method")

    if np.sign(f(a)) == 0:
        return a
    if np.sign(b) == 0:
        return b

    xi = (a + b) / 2

    if np.sign(f(a)) * np.sign(f(xi)) < 0:
        if abs((xi - a) / xi) < e:
            return xi
        else:
            return bisection(f, a, xi, it, e)
    elif np.sign(f(xi)) * np.sign(f(b)) < 0:
        if abs((xi - b) / xi) < e:
            return xi
        else:
            return bisection(f, xi, b, it, e)
    


