import numpy as np


from la.linear_algebra import *
from la.linalg import *



def print_mat(matrix: list, digits: int = 3):
    for row in matrix:
        for col in row:
            print(round(col, digits), end='\t')
        print()


if __name__ == "__main__":

    matrix = [
        [1.0, 2.0, 3.0],
        [4.0, 5.0, 6.0]
    ]

    v = [1, 2, 3]
    # print(f'Norm of the vector: {v} = ', end='')
    # print(norm(v))
    # print(f'Norm of the matrix: ')
    # print_mat(matrix)   
    # print(norm(matrix))

    transpose_matrix = transpose(matrix)

    # print('Initial matrix')
    # print_mat(matrix)
    # print('Transpose matrix')
    # print_mat(transpose_matrix)
    # print()

    m1 = [
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9]
    ]
    m2 = [
        [0.1, 0.2, 0.3],
        [0.4, 0.5, 0.6],
        [0.7, 0.8, 0.9]
    ]
    # print_mat(mat_add(m1, m2))
    # print()
    # print_mat(mat_sub(m1, m2))
    # print()
    # print_mat(scalar_mul(2, m2))
    # print()

    m1 = [
        [1.0, 4.0],
        [2.0, 2.0],
        [3.0, 5.0]
    ]
    m2 = [
        [2.0, 3.0, 1.0, 5.0],
        [0.0, -1.0, 2.0, 4.0]
    ]

    m3 = mat_mul(m1, m2)

    m4 = [
        [1.0, 2.0, 3.0],
        [1.0, 1.5, 2.5],
        [0.1, 0.5, 1.5]
    ]
    m5 = inverse(m4)

    # print_mat(m4)
    # print()
    # print_mat(m5)
    # print('Check inversed matrix')
    # print_mat(mat_mul(m5, m4))
    # print('Check inversed matrix')
    # print_mat(mat_mul(m4, m5))

    print()
    matrix = [
        [3, 2, -3],
        [7, -1, 0],
        [2, -4, 5]
    ]
    print(det(matrix))
    print(np.linalg.det(matrix))

    # Check is the matrix diagonal dominant
    a = [
        [3.0, -2.0, 1.0],
        [1.0, 3.0, 2.0],
        [-1.0, 2.0, 4.0]
    ]
    b = [
        [-2.0, 2.0, 1.0],
        [1.0, 3.0, 2.0],
        [1.0, -2.0, 0.0]
    ]

    print(is_diagonal_dominant(a), is_diagonal_dominant(b))

    a = [[2, 5, 8, 7], [5, 2, 2, 8], [7, 5, 6, 6], [5, 4, 4, 8]]

    mat_pow = np.linalg.matrix_power(a, 5)
    for pi in mat_pow:
        print(pi)

    p, l, u = lup(a)
    for pi in p:
        print(pi)
    print('L')
    for li in l:
        print([round(lij, 4) for lij in li])
    print('U')
    for ui in u:
        print([round(uij, 4) for uij in ui])
    print()
