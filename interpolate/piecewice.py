import os
import sys
import scipy


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

    def __init__(self, coef: list):
        self.__degree = len(coef)
        self.__coef = coef

    @property
    def degree(self) -> int:
        return self.__degree

    def __call__(self, var: float | int) -> float | int:
        return sum([c * var ** i for i, c in enumerate(self.__coef)])


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
        elif self.__n == 2:
            # compute coefficients
            pts = deepcopy(self.__points)
            y = [0.0 for row in range(3 * (k - 1))]
            x = [[0.0 for col in range(3 * (k - 1))] for row in range(3 * (k - 1))]

            xrow = (k - 2) * 2 + 1

            # Setup [X] and [Y]
            y[0] = pts[0][1]
            y[xrow] = pts[-1][1]
            
            x[0][0], x[0][1], x[0][2] = 1, pts[0][0], pts[0][0] ** 2
            x[xrow][xrow + 1], x[xrow][xrow + 2], x[xrow][xrow + 3] = 1, pts[-1][0], pts[-1][0] ** 2
            x[-1][2] = 2

            rows = k - 2
            for row in range(1, rows + 2, 2):
                xrow = row + k
                if row == 1:
                    xi = pts[row][0]
                    yi = pts[row][1]
                    
                    x[row][0], x[row][1], x[row][2] = 1, xi, xi ** 2
                    x[row + 1][row + 2], x[row + 1][row + 3], x[row + 1][row + 4] = 1, xi, xi ** 2
                    x[xrow + 1][1], x[xrow + 1][2] = 1, 2 * xi
                    x[xrow + 1][4], x[xrow + 1][5] = -1, -2 * xi

                    y[row], y[row + 1] = yi, yi
                else:
                    xi = pts[row - 1][0]
                    yi = pts[row - 1][1]
                    x[row][row], x[row][row + 1], x[row][row + 2] = 1, xi, xi ** 2
                    x[row + 1][row + 3], x[row + 1][row + 4], x[row + 1][row + 5] = 1, xi, xi ** 2
                    x[xrow][row], x[xrow][row + 1], x[xrow][row + 2] = 0.0, 1.0, 2 * xi
                    x[xrow][row + 3], x[xrow][row + 4], x[xrow][row + 5] = 0.0, 1.0, 2 * xi
                    y[row], y[row + 1] = yi, yi

            print(x)
            
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
