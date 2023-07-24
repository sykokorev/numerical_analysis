

def norm(v: list) -> float:
    if not isinstance(v, list):
        raise TypeError('Vector must be list type')    
    
    return sum([vi ** 2 for vi in v]) ** 0.5


def distance(v1: list, v2: list) -> float:
    if not isinstance(v1, list) or not isinstance(v2, list):
        raise TypeError('Vectors must be list type')    
    
    if len(v1) != len(v2):
        raise ValueError('Vectors must be same dimention')
    
    return sum([(v1i - v2i) ** 2 for v1i, v2i in zip(v1, v2)]) ** 0.5


def inverse(mat: list) -> list:

    rows, cols = len(mat), len(mat[0])
    res = [[0.0 for i in range(rows)] for j in range(cols)]
    for r in range(rows):
        for c in range(cols):
            res[c][r] = mat[r][c]
    return res


def mat_mul(m1: list, m2: list) -> list:

    n, m, p, q = len(m1), len(m1[0]), len(m2), len(m2[0])
    if m != p:
        raise ValueError('Number of rows of first matrix mus be equal number of columns second matrix')
    
    res = [[0.0 for i in range(q)] for j in range(n)]
    for i in range(n):
        for j in range(q):
            sum = 0.0
            for k in range(p):
                sum += m1[i][k]  * m2[k][j]
            res[i][j] = sum

    return res
