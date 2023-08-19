import numpy as np
import scipy

from la.linear_algebra import *
from la.linalg import *


if __name__ == "__main__":

    A = [
        [-2.0, 2.0, -1.0],
        [6.0, -6.0, 7.0],
        [3.0, -8.0, 4.0]
    ]

    print('Init matrix')
    for ai in A:
        print(ai)

    P, L, U = lup(A, permute_l=False)

    print('my')
    print('P')
    for pi in P:
        print(*pi)

    print('L')
    for li in L:
        print(str([round(lij, 2) for lij in li]))

    print('U')
    for ui in U:
        print(str([round(uij, 2) for uij in ui]))

    print("Checking")
    for mi in mat_mul(mat_mul(P, L), U):
        print(*mi)

    P, L, U = scipy.linalg.lu(A)

    print('scipy')
    print('P')
    for pi in P:
        print(*pi)

    print("L")
    for li in L:
        print(str([round(lij, 2) for lij in li]))

    print('U')
    for ui in U:
        print(str([round(uij) for uij in ui]))


    # X = [
    #         [1, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
    #         [1, -0.8, 0.64, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
    #         [0.0, 0.0, 0.0, 1, -0.8, 0.64, 0.0, 0.0, 0.0], 
    #         [0.0, 0.0, 0.0, 1, -0.6, 0.36, 0.0, 0.0, 0.0], 
    #         [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1, -0.6, 0.36], 
    #         [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1, -0.4, 0.16], 
    #         [0.0, 1, -1.6, 0.0, -1, 1.6, 0.0, 0.0, 0.0], 
    #         [0.0, 0.0, 0.0, 0.0, 1.0, -1.2, 0.0, 1.0, -1.2], 
    #         [0.0, 0.0, 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    #     ]
    

    # P, L, U = lup(X)

    # print('my')
    # for pi in P:
    #     print(*pi)

    # # Checking
    # # X = mat_mul(mat_mul(P, L), U)

    # # for xi in X:
    # #     print(str([round(xij, 2) for xij in xi]))

    # P, L, U = scipy.linalg.lu(X)

    # print('scipy')
    # for pi in P:
    #     print(*pi)
