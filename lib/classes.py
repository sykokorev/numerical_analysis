from __future__ import annotations

import os
import sys
import math

from copy import deepcopy
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

    def  __add__(self, f: Function) -> Addition:
        return Addition([self, f])
    
    def __radd__(self, f: Function) -> Addition:
        return Addition([self, f])
    
    def __iadd__(self, f: Function) -> Addition:
        return Addition([self, f])
    
    def __mul__(self, f: Function) -> Product:
        return Product([self, f])
    
    def __rmul__(self, f) -> Product:
        return Product([self, f])
    
    def __imul__(self, f) -> Product:
        return Product([self, f])


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
class Product:
    func: List[Function]

    def __call__(self, *args: Any, **kwds: Any) -> Any:
        res = 1.0
        for f in self.func:
            res *= f(args[0])
        return res

    def __str__(self):
        s = 'f(x) = ' + str(self.func[0])
        for f in self.func[1:]:
            s += ' * ' + str(f)
        return s
    
    def derivative(self):
        tmp = []
        # for i, f in enumerate(self.func):
        #     tmp.append([f.derivative(), self.func[j] for j in range(len(self.func)-1) if j != i])
        # print(tmp)


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
        a = self.a if self.a != 1.0 else ''
        return f'{a} * x^{self.degree}'


@dataclass
class Sin(Function):
    '''
        f(x) = a * sin(b * x)
    '''
    a: float = 1.0
    b: float = 1.0

    def __call__(self, *args: Any, **kwds: Any) -> Any:
        return self.a * math.sin(args[0] * self.b)

    def derivative(self) -> Cos:
        return Cos(self.a, self.b)

    def differential(self, x: float) -> float:
        c = self.derivative()
        return c(x)

    def __str__(self):
        a = self.a if self.a != 1.0 else ''
        b = self.b if self.b != 1.0 else ''
        return f'{a}sin({b}x)'


@dataclass
class Cos(Function):
    '''
        f(x) = a * cos(b * x)
    '''    
    a: float = 1.0
    b: float = 1.0

    def __call__(self, *args: Any, **kwds: Any) -> Any:
        
        return self.a * math.cos(args[0] * self.b)

    def derivative(self) -> Sin:
        return Sin((-1) * self.a, self.b)

    def differential(self, x: float) -> float:
        s = self.derivative()
        return s(x)

    def __str__(self):
        a = self.a if self.a != 1.0 else ''
        b = self.b if self.b != 1.0 else ''
        return f'{a}cos({b}x)'
