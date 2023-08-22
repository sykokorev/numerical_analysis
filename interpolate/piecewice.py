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

    @property
    def intervalX(self):
        return self.__interval

    @intervalX.setter
    def intervalX(self, x: list):
        if isinstance(x, list) and len(x) == 2:
            if any([isinstance(xi, (float, int)) for xi in x]):
                self.__interval = x

    def eval(self):
        self.__a, self.__b = mat_mul(
            [self.__p1[1], self.__p2[1]],
            inverse([[self.__p1[0], 1], [self.__p2[0], 1]])
        )
        return None

    def __call__(self, x: float):
        return (x, self.__a * x + self.__b)

    def __str__(self):
        return f'Line: [p1={self.__p1}, p2={self.__p2}, a={self.__a}, b={self.__b}]'


class Polynomial:

    def __init__(self, coef: list):
        self.__degree = len(coef)
        self.__coef = coef

    @property
    def degree(self) -> int:
        return self.__degree - 1

    @property
    def coefficients(self):
        return self.__coef

    def __call__(self, x: float | int) -> float | int:
        return (x ,sum([c * x ** i for i, c in enumerate(self.__coef)]))

    def __str__(self):
        return f'Polynomial degree {self.__degree - 1}, coefficients: {self.__coef}'

class Interpolate:

    def __init__(self, points: list, n: int, spline: str='not_a_knot'):

        self.__n = n
        self.__objects = None
        self.__points = points
        self.__intervals = None
        self.__spline = spline
    
    @property
    def n(self):
        return self.__n

    @n.setter
    def n(self, n: int):

        if not isinstance(n, int): return None
        self.__n = n
        self.eval()

    @property
    def intervals(self):
        return self.__intervals

    @property
    def functions(self):
        return self.__objects

    def __line(self):
        self.__objects: List[Line] = []
        if self.__n == 1:
            for i in range(k - 1):
                self.__objects.append(Line(
                    self.__points[i], self.__points[i+1]
                ))
                self.__objects[i].eval()
        self.__intervals = []
        for pt in self.__points:
            self.__intervals.append(pt[0])

    def __quadratic(self):
        k = len(self.__points)
        self.__objects: List[Polynomial] = []
        pts = deepcopy(self.__points)
        total = (self.__n + 1) * (k - 2)
        x = zeros(total, axis=2)
        y = zeros(total)
        
        # Equations rows
        rows = k
        for j in range(self.__n + 1):
            x[0][j] = pts[0][0] ** j
            y[0] = pts[0][1]

        col = 0
        for i, row in enumerate(range(1, k - 1), 1):
            xi = pts[i][0]
            yi = pts[i][1]
            for j in range(self.__n + 1):
                x[row][col + j] = xi ** j
            y[row] = yi
            col += self.__n + 1

        col = total - 1 - self.__n
        for j in range(self.__n + 1):
            x[rows - 1][col + j] = pts[-1][0] ** j
        y[rows - 1] = pts[-1][1]
        # First order derivatives rows. To interval (xi + x(i+1) / 2)
        frow = rows
        srow = k - 3
        rows = 2 * (k - 3)
        col = 0
        for i, row in enumerate(range(frow, frow + k - 3), 1):
            x1, x2 = pts[i][0], pts[i + 1][0]
            xi = (x1 + x2) / 2
            for j in range(self.__n + 1):
                x[row][col + j] = xi ** j
                x[row][col + j + self.__n + 1] = (-1) * xi ** j
                if xi:
                    x[srow + row][col + j] = j * xi ** (j - 1)
                    x[srow + row][col + j + self.__n + 1] = (-1) * j * xi ** (j - 1)
                elif xi == 0.0:
                    x[srow + row][col + j] = 1.0 if j == 1 else 0.0
                    x[srow + row][col + j + self.__n + 1] = -1.0 if j == 1 else 0.0

            col += self.__n + 1
        
        coef = self.__comp_coef(x, y)
        for i in range(0, len(coef) - 1, self.__n + 1):
            self.__objects.append(
                Polynomial(coef=coef[i:i + self.__n + 1])
            )
        self.__intervals = []
        for pt in self.__points:
            self.__intervals.append(pt[0])

    def __cubic(self):
        
        k = len(self.__points)
        self.__objects: List[Polynomial] = []
        pts = deepcopy(self.__points)

        if self.__spline == 'natural':
            total = (self.__n + 1) * (k - 1)
            x = zeros(total, axis=2)
            y = zeros(total)

            # Equation rows si-1(xi) = si(xi) = yi
            rows = 2 * k - 4
            for j in range(self.__n + 1):
                x[0][j] = pts[0][0] ** j
            y[0] = pts[0][1]

            col = self.__n + 1
            for i, row in enumerate(range(1, rows, 2), 1):
                xi = pts[i][0]
                yi = pts[i][1]
                for j in range(self.__n + 1):
                    x[row][j + col - (self.__n + 1)] = xi ** j
                    x[row + 1][j + col] = xi ** j
                y[row], y[row + 1] = yi, yi
                col += self.__n + 1
            
            col = total - (self.__n + 1)
            for j in range(self.__n + 1):
                x[rows + 1][j + col] = pts[-1][0] ** j
            y[rows + 1] = pts[-1][1]

            col = 0
            frow = rows + 2
            rows = frow + k - 2
            # First derivative rows
            for i, row in enumerate(range(frow, rows), 1):
                xi = pts[i][0]
                for j in range(self.__n + 1):
                    if xi:
                        x[row][col + j] = j * xi ** (j - 1)
                        x[row][col + j + self.__n + 1] = (-1) * j * xi ** (j - 1)
                    elif xi == 0.0:
                        x[row][col + j] = 1.0 if j == 1 else 0.0
                        x[row][col + j + self.__n + 1] = -1.0 if j == 1 else 0.0
                col += self.__n + 1
            
            col = 0
            frow = rows
            rows = frow + k - 2
            # Second derivative rows
            for i, row in enumerate(range(frow, rows), 1):
                xi = pts[i][0]
                for j in range(self.__n + 1):
                    if xi:
                        x[row][j + col] = j * (j - 1) * xi ** (j - 2)
                        x[row][j + col + self.__n + 1] = (-1) * j * (j - 1) * xi ** (j - 2)
                    elif xi == 0.0:
                        x[row][col + j] = 1.0 if j == 2 else 0.0
                        x[row][col + j + self.__n + 1] = -1.0 if j == 2 else 0.0
                col += self.__n + 1

                x[-2][2], x[-2][3] = 2, 6 * pts[0][0]
                x[-1][-2], x[-1][-1] = 2, 6 * pts[-1][0]

            coef = self.__comp_coef(x, y)
            for i in range(0, len(coef) - 1, self.__n + 1):
                self.__objects.append(Polynomial(coef=coef[i : i + self.__n + 1]))
            
            self.__intervals = []
            for pt in self.__points:
                self.__intervals.append(pt[0])

        elif self.__spline == 'not_a_knot':
            total = (self.__n + 1) * (k - 3)
            x = zeros(total, axis=2)
            y = zeros(total)

            # First two rows
            for j in range(self.__n + 1):
                x[0][j] = pts[0][0] ** j
                x[1][j] = pts[1][0] ** j
            y[0] = pts[0][1]
            y[1] = pts[1][1]

            # Equations rows
            col = 0
            frow = 2
            rows =  2 * (k - 3)
            for i, row in enumerate(range(frow, rows, 2), 2):
                xi, yi = pts[i][0], pts[i][1]
                for j in range(self.__n + 1):
                    x[row][j + col] = xi ** j
                    x[row + 1][j + col + self.__n + 1] = xi ** j
                y[row], y[row + 1] = yi, yi
                col += self.__n + 1

            frow = 2 * (k - 3)
            for j in range(self.__n + 1):
                x[rows][col + j] = pts[-2][0] ** j
                x[rows + 1][col + j] = pts[-1][0]  ** j
            y[rows] = pts[-2][1]
            y[rows + 1] = pts[-1][1]

            # First derivatives rows
            col = 0
            frow = 2 * (k - 2)
            rows = frow + (k - 4)
            for i, row in enumerate(range(frow, rows), 2):
                xi = pts[i][0]
                for j in range(self.__n + 1):
                    if xi:
                        x[row][col + j] = j * xi ** (j - 1)
                        x[row][col + j + self.__n + 1] = (-1) * j * xi ** (j - 1)
                    else:
                        x[row][col + j] = 1.0 if j == 2 else 0.0
                        x[row][col + j + self.__n + 1] = -1.0 if j == 2 else 0.0
                col += self.__n + 1

            # Second derivatives rows
            col = 0
            frow = 2 * (k - 2) + (k - 4)
            rows = frow + (k - 4)
            for i, row in enumerate(range(frow, rows), 2):
                xi = pts[i][0]
                for j in range(self.__n + 1):
                    if xi:
                        x[row][col + j] = j * (j - 1) * xi ** (j - 2)
                        x[row][col + j + self.__n + 1] = (-1) * j * (j - 1) * xi ** (j - 2)
                    else:
                        x[row][col + j] = 1.0 if j == 2 else 0.0
                        x[row][col + j + self.__n + 1] = -1.0 if j == 2 else 0.0
                col += self.__n + 1

            coef = self.__comp_coef(x, y)
            for i in range(0, len(coef) - 1, self.__n + 1):
                self.__objects.append(Polynomial(coef=coef[i : i + self.__n + 1]))

            self.__intervals = []
            for i, pt in enumerate(self.__points):
                if (i != 1 and i != len(self.__points) - 2):
                    self.__intervals.append(pt[0])

    def __comp_coef(self, A, b):
        return solveLUP(A, b)
    
    def eval(self, n=None, spline='not_a_knot'):
        if n and isinstance(n, float | int): self.__n = n

        k = len(self.__points)
        if self.__n == 1:
            self.__line()
        elif self.__n == 2:
            self.__quadratic()
        elif self.__n == 3:
            self.__cubic()

    def __call__(self, x):
        idx = self.__index(x)
        return self.__objects[idx](x)

    def __index(self, x):
        xs = [xi for xi in self.__intervals]
        n = len(xs)
        for i in range(n-1):
            if x >= xs[i] and x <= xs[i+1]:
                return i
        return i

    def __str__(self):
        out = ''
        for i, obj in enumerate(self.__objects, 1):
            out += f'{str(obj)}, interval {self.__intervals[i-1:i+1]}\n'
        return out
