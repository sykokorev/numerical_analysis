import sys


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
