from __future__ import annotations

from typing import Any
import numpy as np
import matplotlib.pyplot as plt

from dataclasses import dataclass
from abc import ABCMeta, abstractmethod


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


if __name__ == "__main__":

    func = Exponential(0.5, 1)
    
    print(taylor_series(func, 5, 4, 1))
