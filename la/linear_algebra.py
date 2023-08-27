from copy import deepcopy


# Swap row i and j
def swap_rows(A, i, j):

    a = deepcopy(A)
    if i == j: return a
    tmp = [ai for ai in a[j]]
    a[j] = [ai for ai in a[i]]
    a[i] = [aj for aj in tmp]
    return a


def swap_cols(A, i, j):

    a = deepcopy(A)
    rows = len(a)
    if i == j: return a
    tmp = []
    for row in range(rows):
        tmp.append(a[row][j])
    for row in range(rows):
        a[row][j] = a[row][i]
    for row in range(rows):
        a[row][i] = tmp[row]
    return a


def max_element(A, row, col):
    n = len(A)
    idx = (row, col)
    max_arg = abs(A[row][col])
    for r in range(row, n):
        for c in range(col, n):
            if (narg:=abs(A[r][c])) >= max_arg:
                max_arg = narg
                idx = (r, c)
    
    return idx, max_arg


def norm(v: list) -> float:
    if not isinstance(v, list):
        raise TypeError('Vector must be list type')    
    
    if all([isinstance(vi, (int, float)) for vi in v]):
        return sum([vi ** 2 for vi in v]) ** 0.5
    elif all([isinstance(vi, list) for vi in v]):
        s = 0.0
        for row in v:
            s += sum([col ** 2 for col in row])
        return s ** 0.5


def distance(v1: list, v2: list) -> float:
    if not isinstance(v1, list) or not isinstance(v2, list):
        raise TypeError('Vectors must be list type')    
    
    if len(v1) != len(v2):
        raise ValueError('Vectors must be same dimention')
    
    return sum([(v1i - v2i) ** 2 for v1i, v2i in zip(v1, v2)]) ** 0.5


def dot(m1: list, m2: list) -> float or bool:

    if len(m1) != len(m2):
        raise ValueError('Shapes of matrices not aligned')
    
    try:
        if len(m1[0]) != len(m2[0]):
            return ValueError('Shape of matrices not aligned')
        elif len(m1[0]) != len(m1):
            raise ValueError('Shapes of fmatrices not aligned')
        col = len(m1[0])        
    except TypeError:
        col = 0

    if not col:
        return sum([m1i * m2i for (m1i, m2i) in zip(m1, m2)])
    else:
        out = 0.0
        for m1i, m2i in zip(m1, m2):
            out += sum(m1ij * m2ij for (m1ij, m2ij) in zip(m1i, m2i))
        return out


def mat_mul(m1: list, m2: list) -> list:

    # Multiply a matrix by a scalar
    if isinstance(m1, int | float) and hasattr(m2, '__iter__'):
        if isinstance(m2[0], int | float): 
            return [m1 * m2i for m2i in m2]
        else:
            return [mat_mul(m1, m2[i]) for i in range(len(m2))]

    if isinstance(m2, int | float) and hasattr(m1, '__iter__'):
        if isinstance(m1[0], int | float): 
            return [m2 * m1i for m1i in m1]
        else:
            return [mat_mul(m2, m1[i]) for i in range(len(m1))]

    # Multiply a matrix by a vector

    if isinstance(m1[0], (float, int)) and hasattr(m2[0], '__iter__'):
        res = [0.0 for i in range(len(m2))]
        for i, el in enumerate(m2):
            res[i] = sum([eli * m1i for eli, m1i in zip(el, m1)])
        return res
        
    if isinstance(m2[0], (float, int)) and hasattr(m1[0], '__iter__'):
        res = [0.0 for i in range(len(m1))]
        for i, el in enumerate(m1):
            res[i] = sum([eli * m2i for eli, m2i in zip(el, m2)])
        return res

    # Multiply a matrix by a mtrix

    n, m, p, q = len(m1), len(m1[0]), len(m2), len(m2[0])

    if m != p:
        raise ValueError('Number of rows of first matrix mus be equal number of columns second matrix')
    
    res = [[0.0 for i in range(q)] for j in range(n)]
    for i in range(n):
        for j in range(q):
            summ = 0.0
            for k in range(p):
                summ += m1[i][k]  * m2[k][j]
            res[i][j] = summ

    return res


def scalar_mul(m1: list | float, m2: list | float) -> list:

    if isinstance(m1, list) and isinstance(m2, (int, float)):
        if isinstance(m1[0], (float, int)):
            return [m2 * mi for mi in m1]
        elif isinstance(m1[0], list):
            return [scalar_mul(mi, m2) for mi in m1]
    
    if isinstance(m2, list) and isinstance(m1, (int, float)):
        if isinstance(m2[0], (float, int)):
            return [m1 * mi for mi in m2]
        elif isinstance(m2[0], list):
            return [scalar_mul(mi, m1) for mi in m2]


def mat_add(m1: list, m2: list) -> list:

    if isinstance(m1, int | float) and isinstance(m2, int | float):
        return m1 + m2

    if isinstance(m1, int | float) and hasattr(m2, '__iter__'):
        if isinstance(m2[0], int | float): 
            return [m1 + m2i for m2i in m2]
        else:
            return [mat_add(m1, m2[i]) for i in range(len(m2))]

    if isinstance(m2, int | float) and hasattr(m1, '__iter__'):
        if isinstance(m1[0], int | float): 
            return [m2 + m1i for m1i in m1]
        else:
            return [mat_add(m2, m1[i]) for i in range(len(m1))]
    
    if hasattr(m1, '__iter__') and hasattr(m2, '__iter__'):
        if isinstance(m1[0], int | float) and isinstance(m2[0], int | float):
            return [m1i + m2i for m1i, m2i in zip(m1, m2)]
        elif hasattr(m1[0], '__iter__') and hasattr(m2[0], '__iter__'):
            return [mat_add(m1[i], m2[i]) for i in range(len(m1))]


def mat_sub(m1: list, m2: list) -> list:
    if isinstance(m1, int | float) and hasattr(m2, '__iter__'):
        if isinstance(m2[0], int | float): 
            return [m1 - m2i for m2i in m2]
        else:
            return [mat_sub(m1, m2[i]) for i in range(len(m2))]

    if isinstance(m2, int | float) and hasattr(m1, '__iter__'):
        if isinstance(m1[0], int | float): 
            return [m2 - m1i for m1i in m1]
        else:
            return [mat_sub(m2, m1[i]) for i in range(len(m1))]


def transpose(m) -> list:
    return list(map(list, zip(*m)))


def separate(m: list, col_num: int) -> list:
    nrows = len(m)
    cols = len(m[0])
    ncols1 = col_num
    ncols2 = cols - col_num

    m1 = [[0] * ncols1 for row in range(nrows)]
    m2 = [[0] * ncols2 for row in range(nrows)]

    for row in range(nrows):
        for col in range(cols):
            if (col < col_num):
                m1[row][col] = m[row][col]
            else:
                m2[row][col - col_num] = m[row][col]

    return [m1, m2]


def is_square(m) -> bool:
    if len(m) == 1 and isinstance(m[0], (int, float)):
        return True
    if not (len(m) == len(m[0])):
        return False
    else:
        return True


def is_diagonal_dominant(m: list) -> bool:

    n = len(m)
    r = []
    for i in range(n):
        r.append(abs(m[i][i]) >= sum([m[i][j] for j in range(n) if j != i]))
    return all(r)


def swap_row(m: list, i: int, j: int) -> list:
    temp = [el for el in m[i]]
    m[i] = m[j]
    m[j] = temp
    return m


def join(m1, m2) -> list:
    out = []
    for row1, row2 in zip(m1, m2):
        tmp = [r for r in row1]
        for val in row2:
            tmp.append(val)
        out.append(tmp)

    return out



def upper(m):

    rows = len(m)
    c = deepcopy(m)
    cols = len(c[0])
    
    for k in range(rows):
        pivot = c[k][k]
        for i in range(k+1, rows):
            pivot_col = c[i][k]
            for j in range(cols):
                pivot_row = c[k][j]
                c[i][j] = c[i][j] - pivot_col * pivot_row / c[k][k]

    return c


def joint(a, b):

    out = []
    rows = len(a)
    for row in range(rows):
        out.append([ai for ai in a[row]] + [b[row]])
    
    return out


def split(c: list, col: int = 1):
    a = []
    b = []
    for row in c:
        a.append(row[:-col])
        tmp = row[col+2:]
        b.append(tmp[0] if len(tmp) == 1 else tmp)
    return (a, b)


def find_row_with_max_element(m, col_num: int=0, starting_row: int=0):
    tmp = m[starting_row][col_num]
    row_idx = starting_row
    for k in range(starting_row+1, len(m)):
        if abs(m[k][col_num]) > abs(tmp):
            row_idx = k
            tmp = m[k][col_num]
    return row_idx


def identity(ndim: int) -> list:
    
    return [[float(i == j) for i in range(ndim)] for j in range(ndim)]


def zeros(ndim: int, axis: int = 1) -> list:
    if axis == 1:
        return [0.0 for i in range(ndim)]
    elif axis == 2:
        return [[0.0 for i in range(ndim)] for j in range(ndim)]
    else:
        return None


def compare(m1: list, m2: list) -> bool:
    
    if not (len(m1) - len(m2)) or not (len(m1[0]) - len(m2[0])):
        return False
        
    for m1_row, m2_row in zip(m1, m2):
        for el_1, el_2 in zip(m1_row, m2_row):
            if (el_1 - el_2) > 10 ** -6:
                return False

    return True


def inverse(m, diag: bool = False) -> list:
    '''
        m - Matrix must be square
        diag - is the m diagonal matrix
    '''

    if not is_square(m):
        return False

    identity_matrix = identity(ndim=len(m))
    num_rows = len(m) - 1
    joint_matrix = join(m, identity_matrix)

    flag = False
    count = 0
    max_count = 100

    if diag:
        inverse_matrix = deepcopy(m)
        if all([inverse_matrix[i][i] != 0 for i in range(len(inverse_matrix))]):
            for i in range(len(inverse_matrix)):
                inverse_matrix[i][i] = 1 / inverse_matrix[i][i]
            return inverse_matrix
        else:
            raise ZeroDivisionError('All main diagonal elements must be non-zero.')

    while not flag and count < max_count:
        for i in range(num_rows + 1):
            if joint_matrix[i][i] == 0.0:
                max_el_idx = find_row_with_max_element(joint_matrix, i, i)
                joint_matrix = swap_row(joint_matrix, i, max_el_idx)
            div_e = joint_matrix[i][i]
            factor = 1 / div_e
            joint_matrix[i] = [e * factor for e in joint_matrix[i]]
            for row in range(0, num_rows+1):
                if row != i:
                    if joint_matrix[row][i] != 0:
                        sub = (-1) * joint_matrix[row][i]
                        row_add = [el * sub for el in joint_matrix[i]]
                        joint_matrix[row] = [e1 + e2 for e1, e2 in zip(row_add, joint_matrix[row])]
    
        identity_matrix, inverse_matrix = separate(m=joint_matrix, col_num=num_rows+1)
        if not compare(identity(ndim=num_rows+1), identity_matrix):
            flag = True
        count += 1
    
    if not flag:
        return False
    else:
        return inverse_matrix


def det(m: list, mul=1.0) -> float | bool:

    if not is_square(m):
        return False

    width = len(m)

    if width == 1:
        return mul * m[0][0]
    else:
        sign = -1
        answer = 0
        for col in range(width):
            mi = []
            for row in range(1, width):
                buff = []
                for k in range(width):
                    if k != col:
                        buff.append(m[row][k])
                mi.append(buff)
            sign *= -1
            answer = answer + mul * det(mi, sign * m[0][col])
    return answer

