import sys

from typing import Tuple, List


from la.linear_algebra import *


def cramer(a, b) -> list | bool:
    
    if not is_square(a):
        return False
    
    detA = det(a)
    if not detA: 
        return None

    x = []
    width = len(a)

    for col in range(width):
        d = deepcopy(a)

        for row in range(width):
            d[row][col] = b[row]

        x.append(det(d) / detA)
    return x


def naive_gaussian(a, b):

    ai, bi = split(upper(joint(a, b)))
    rows = len(ai)
    x = [0.0 for i in range(rows)]

    for i in range(rows-1, -1, -1):
        x[i] = (1 / ai[i][i]) * (bi[i] - sum([ai[i][j] * x[j] for j in range(i + 1, rows)]))
    return x


def gauss_jordan(a, b):

    c = joint(a, b)
    rows = len(c)
    cols = len(c[0])

    for i in range(rows):

        if c[i][i] == 0.0:
            sys.exit('Divide by zero detected')

        for j in range(rows):
            if i != j:
                ratio = c[j][i] / c[i][i]

                for k in range(rows + 1):
                    c[j][k] = c[j][k] - ratio * c[i][k]

    for k in range(rows):
        pivot = c[k][k]
        for j in range(cols):
            c[k][j] /= pivot

    x = [ci[-1] for ci in c]

    return x


# LU Decomposition. Return tuple(L, U)
def lu(a):

    width = len(a)
    upper = deepcopy(a)
    lower = identity(width)

    for k in range(width):
        
        pivot = upper[k][k]
        for i in range(k+1, width):
            if k != i:
                ratio = upper[i][k] / pivot
                for j in range(width):
                    upper[i][j] = upper[i][j] - upper[k][j] * ratio
                lower[i][k] = ratio

    return (lower, upper)    


def solveLU(a, b):

    width = len(a)
    l, u = lu(a)
    y = [0.0 for i in range(width)]
    x = [0.0 for i in range(width)]

    # Forward substitution [L]{y} = {b}, {y} = [U]{x}
    for i in range(width):
        y[i] = (1 / l[i][i]) * (b[i] - sum(l[i][j] * y[j] for j in range(i)))

    # Backward substitution [L]{x} = {y}
    for i in range(width-1, -1, -1):
        x[i] = (1 / u[i][i]) * (y[i] - sum(u[i][j] * x[j] for j in range(i+1, width)))

    return x


def cholesky(a):

    width = len(a)
    u = [[0.0 for i in range(width)] for j in range(width)]

    for i in range(width):
        u[i][i] = (a[i][i] - sum([u[k][i] ** 2 for k in range(i)])) ** 0.5
        for j in range(i+1, width):
            u[i][j] = (1 / u[i][i]) * (a[i][j] - sum([u[k][i] * u[k][j] for k in range(i)]))
   
    return u


def cholab(a, b):

    width = len(a)
    u = cholesky(a)
    ut = transpose(u)
    y = [0.0 for i in range(width)]
    x = [0.0 for i in range(width)]

    # Forward substitution 
    for i in range(width):
        y[i] = (1 / ut[i][i]) * (b[i] - sum(ut[i][j] * y[j] for j in range(i)))

    # Backward substitution
    for i in range(width-1, -1, -1):
        x[i] = (1 / u[i][i]) * (y[i] - sum(u[i][j] * x[j] for j in range(i+1, width)))
    
    return x