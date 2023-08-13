import os
import sys


from copy import deepcopy


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(DIR)


from la.linalg import *
from la.linear_algebra import *


class Line:

    def __init__(self, p1, p2):
        self.__p1 = deepcopy(p1)
        self.__p2 = deepcopy(p2)
        self.__a = None
        self.__b = None

    @property
    def p1(self):
        return self.__p1

    @property
    def p2(self):
        return self.__p2

    @property
    def a(self):
        return self.__a
    
    @property
    def b(self):
        return self.__b

    def eval(self):
        self.__a, self.__b = mat_mul(
            [self.__p1[1], self.__p2[1]],
            inverse([[self.__p1[0], 1], [self.__p2[0], 1]])
        )
        return None

    def __call__(self, x: float):
        return [x, self.__a * x + self.__b]

    def __str__(self):
        return f'Line: [p1={self.__p1}, p2={self.__p2}, a={self.__a}, b={self.__b}]'


class Polynomial:

    def __init__(self):
        pass


class Interpolate:

    def __init__(self, points: list, n: int):

        self.__n = n
        self.__objects = None
        self.__points = points
    
    @property
    def n(self):
        return self.__n

    @n.setter
    def n(self, n: int):

        if not isinstance(n, int): return None
        self.__n = n
        self.eval()

    def eval(self):
        self.__objects = []
        k = len(self.__points)
        if self.__n == 1:
            for i in range(k - 1):
                self.__objects.append(Line(
                    self.__points[i], self.__points[i+1]
                ))
                self.__objects[i].eval()
        else:
            x = [xi[0] for xi in self.__points]
            print(x)
            y = [yi[1] for yi in self.__points]
            ys = [0.0 for i in range(3*(k-1))]
            xs = [[0.0 for i in range(3 * (k - 1))] for j in range(3 * (k - 1))]
            
            ys[0] = y[0]
            ys[2 * k - 3] = y[-1]

            for i in range(1, k-1):
                for j in range(i, i+2):
                    ys[i + j - 1] = y[i]

            tmp = []
            for i in range(len(x)):
                tmp.append([x[i] ** j for j in range(self.__n + 1)])
            print(tmp)
            for r in range(0, 2 * k - 3, 2):
                xs = [x[i] ** j for j in range(self.__n)]
                

            print(xs)

            
    def __call__(self, x):
        idx = self.__index(x)
        return self.__objects[idx](x)

    def __index(self, x):
        xs = [xi[0] for xi in self.__points]
        n = len(xs)
        for i in range(n-1):
            if x >= xs[i] and x <= xs[i+1]:
                return i
        return i
