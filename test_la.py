from la.linear_algebra import *


def print_mat(matrix: list):
    for row in matrix:
        for col in row:
            print(col, end='\t')
        print()


if __name__ == "__main__":

    matrix = [
        [1.0, 2.0, 3.0],
        [4.0, 5.0, 6.0]
    ]
    
    inverse_matrix = inverse(matrix)

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

    m1_inverse = inverse(m1)
    identity = mat_mul(m1_inverse, m1)
    print_mat(identity)
