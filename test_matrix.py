import numpy as np
import scipy

from scipy.linalg import lu


# from la.linear_algebra import *
# from la.linalg import *


if __name__ == "__main__":
    
    X = [
            [1, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
            [1, -0.8, 0.64, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
            [0.0, 0.0, 0.0, 1, -0.8, 0.64, 0.0, 0.0, 0.0], 
            [0.0, 0.0, 0.0, 1, -0.6, 0.36, 0.0, 0.0, 0.0], 
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1, -0.6, 0.36], 
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1, -0.4, 0.16], 
            [0.0, 1, -1.6, 0.0, -1, 1.6, 0.0, 0.0, 0.0], 
            [0.0, 0.0, 0.0, 0.0, 1.0, -1.2, 0.0, 1.0, -1.2], 
            [0.0, 0.0, 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        ]
    

    P, L, U = lup(X)

    # Checking
    PA = mat_mul(P, X)
    LU = mat_mul(L, U)

    print('[PA]')
    for pai in PA:
        print(pai)

    print('[LU]')
    for lui in LU:
        print(LU)

    P, L, U = lu(X)
    for pi in P:
        print(pi)
