from __future__ import annotations

import os
import sys
import math

from dataclasses import dataclass
from abc import ABCMeta, abstractmethod
from typing import Any, List


LIBDIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__),
    os.pardir)
)
sys.path.append(LIBDIR)


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
class Addition:
    func: List[Function]

    def __call__(self, *args, **kwds: Any) -> Any:
        return sum([f(args[0]) for f in self.func])

    def derivative(self) -> Addition:
        return Addition([f.derivative() for f in self.func])

    def differential(self, *args) -> float:
        return sum([f.derivative()(args[0]) for f in self.func])

    def __str__(self):
        s = 'f(x) = ' + str(self.func[0])
        for f in self.func[1:]:
            s += ' + ' + str(f)
        return s

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

    def __str__(self):
        s = str(self.coeficients[0])
        for i, c in enumerate(self.coeficients[1:], 1):
            if c != 0:
                s += f' + {c}x^{i}'
        return s


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

    def __str__(self):
        return f'{self.a}x^{self.degree}'


@dataclass
class Sin(Function):
    coef: float = 1.0

    def __call__(self, *args: Any, **kwds: Any) -> Any:
        return math.sin(args[0] * self.coef)

    def __add__(self):
        return 

    def derivative(self) -> Cos:
        return Cos(self.coef)

    def differential(self, x: float) -> float:
        c = self.derivative()
        return c(x)

    def __str__(self):
        return f'sin({self.coef}x)'


@dataclass
class Cos(Function):
    coef: float = 1.0

    def __call__(self, *args: Any, **kwds: Any) -> Any:
        
        return math.cos(args[0] * self.coef)

    def derivative(self) -> Sin:
        return Sin(self.coef)

    def differential(self, x: float) -> float:
        s = self.derivative()
        return s(x)

    def __str__(self):
        return f'cos({self.coef}x)'
