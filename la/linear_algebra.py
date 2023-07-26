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


def scalar_mul(m1: list | float, m2: list | float) -> list:

    if isinstance(m1, list) and isinstance(m2, (int, float)):
        return [[m2 * mi for mi in m] for m in m1]
    if isinstance(m2, list) and isinstance(m1, (float, int)):
        return [[m1 * mi for mi in m] for m in m2]
    return None


def mat_add(m1: list, m2: list) -> list:
    return list(map(sum, zip(*i)) for i in zip(m1, m2))


def mat_sub(m1: list, m2: list) -> list:
    return list(map(sum, zip(i, [(-1) * ji for ji in j])) for i, j in zip(m1, m2))


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
    if not (len(m) == len(m[0])):
        return False
    else:
        return True


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


def find_row_with_max_element(m, col_num: int=0, starting_row: int=0):
    tmp = m[starting_row][col_num]
    row_idx = starting_row
    for k in range(starting_row+1, len(m)):
        if abs(m[k][col_num]) > abs(tmp):
            row_idx = k
            tmp = m[k][col_num]
    return row_idx


def identity(ndim: int) -> list:
    
    matrix = [[0.0 for i in range(ndim)] for row in range(ndim)]

    for i in range(len(matrix)):
        matrix[i] = [1 if i == j else 0 for j in range(len(matrix))]
    
    return matrix


def compare(m1: list, m2: list) -> bool:
    
    if not (len(m1) - len(m2)) or not (len(m1[0]) - len(m2[0])):
        return False
        
    for m1_row, m2_row in zip(m1, m2):
        for el_1, el_2 in zip(m1_row, m2_row):
            if (el_1 - el_2) > 10 ** -6:
                return False

    return True


def inverse(m) -> list:
    if not is_square(m):
        return False

    identity_matrix = identity(ndim=len(m))
    num_rows = len(m) - 1
    joint_matrix = join(m, identity_matrix)

    flag = False
    count = 0
    max_count = 100

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
