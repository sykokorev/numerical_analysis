from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import math
import random


from typing import Any
from dataclasses import dataclass
from abc import ABCMeta, abstractmethod


COLORS = [
    'b', 'g', 'r', 'c', 'm', 'y', 'k'
]
MARKERS = [
    '8', '>', '<', '^', 'v', 's',
    'd', 'D', 'H', 'h', '*', 'p'
]


class Function(metaclass=ABCMeta):

    def __init__(self) -> None:
        pass

    @abstractmethod
    def derivative(self) -> Function:
        pass

    @abstractmethod
    def differential(self, x: float) -> float:
        pass


@dataclass
class Polynomial(Function):
    '''
        coefficients: list [a0, a1, a2, ... , an]
        y = a0 + a1 * x + a2 * x ** 2 + a3 * x ** 3 + ... + an * x ** n
    '''
    coeficients: list

    @property
    def degrees(self):
        return len(self.coeficients) - 1
    
    def __call__(self, *args: Any, **kwds: Any) -> Any:
        return sum([ci * args[0] ** i for i, ci in enumerate(self.coeficients)])

    def derivative(self) -> Polynomial:
        return Polynomial([ci * (i + 1) for i, ci in enumerate(self.coeficients[1:])])
    
    def differential(self, x: float) -> float:
        poly = self.derivative()
        return poly(x)


@dataclass
class Exponential(Function):
    '''
        y = a * x ** degree
    '''
    degree: float
    a: float

    def __call__(self, *args: Any, **kwds: Any) -> Any:
        return self.a * args[0] ** self.degree

    def derivative(self) -> Exponential:
        return Exponential(self.degree - 1, self.a * self.degree)
    
    def differential(self, x: float) -> float:
        exp = self.derivative()
        return exp(x)


def fact(n: int):
    
    if n <= 1:
        return 1
    return n * fact(n-1)


def differential(f: Function, x: float, n: int):

    if n <= 0:
        return f(x)

    return differential(f.derivative(), x, n - 1)


def taylor_series(f: function, x: float, x0: float, n: int):

    if n == 0:
        return differential(f, x0, n)
   
    return (differential(f, x0, n) * (x - x0) ** n) / fact(n) + taylor_series(f, x, x0, n - 1)


def max_error(f: Function, x: float, x0: float, n: int) -> list:

    return (max(
        [abs(differential(f, x, n + 1)), 
         abs(differential(f, x0, n + 1))]) * (x - x0) ** (n + 1)) / fact(n + 1)


if __name__ == "__main__":

    # func = Exponential(0.5, 1)
    # n = np.linspace(1, 50, 50)
    
    # errors = []
    # y = []
    # x0 = 8
    # x1 = 5

    # for i in n:
    #     y5 = taylor_series(func, x1, x0, i)
    #     y.append(y5)
    #     errors.append(max_error(func, x1, x0, i))

    # fig, axs = plt.subplots(2, 1)
    # axs[0].plot(n, y, linestyle='-', marker='*', color='red')
    # axs[0].set_title('Function approxiamtion')
    # axs[1].plot(n, errors, linestyle='-', marker='o', color='blue')
    # axs[1].set_title('Residuals')
    # axs[0].grid(True)
    # axs[1].grid(True)

    # plt.show()

    poly = Polynomial([1.2, -0.25, -0.5, -0.15, -0.1])

    fig, ax = plt.subplots()

    x = np.linspace(-2, 2, 21)
    x0 = 0
    n = np.linspace(0, 4, 5)
    y = []

    for xi in x:
        y.append(poly(xi))
    ax.plot(x, y, label='exact', color='k', marker='o')

    for i in n:
        y = []
        for xi in x:
            y.append(taylor_series(poly, xi, x0, i))
        color = random.randint(0, len(COLORS)-1)
        marker = random.randint(0, len(MARKERS) - 1)
        ax.plot(
            x, y, color=COLORS[color], 
            marker=MARKERS[marker], label=f'{int(i)}th order'
        )
        ax

    ax.grid()
    ax.set_xlabel('x')
    ax.set_ylabel('Approximation to f(x) = -0.1$x^4$ - 0.15$x^3$ - 0.5$x^2$ - 0.25x + 1.2')
    fig.legend()
    plt.show()    
