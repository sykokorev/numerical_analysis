import sys

from typing import Tuple, List
from inspect import signature


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

    '''
        |a11  0  ...   0|   | 0  a12 ... a1n|
        |a21 a22 ...   0|   | 0   0  ... a2n|
        |...............|   | ............. |
        |an1 an2 ...   0|   | 0   0  ...  0 |
    '''

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


def pivot_mat(A):

    n = len(A)

    id_mat = identity(n)

    for j in range(n):
        row = j
        max_arg = abs(A[j][j])
        for r in range(j, n):
            if abs(A[r][j]) >= max_arg:
                row = r

        if j != row:
            id_mat[j], id_mat[row] = id_mat[row], id_mat[j]
        
    return id_mat


def lup(A, permute_l = False):

    n = len(A) # Matrix width
    P = identity(n) # Pivot matix. Storage permutations.
    U = deepcopy(A) # Upper triangular matrix
    L = identity(n) # lower triangular matrix

    # Partial pivoting
    for k in range(n):
        j = max(range(k, n), key=lambda i: abs(U[i][k]))
        I = identity(n)
        I[k], I[j] = I[j], I[k]
        P = mat_mul(P, I)
        U[k], U[j] = U[j], U[k]

        for row in range(k + 1, n):
            if U[row][k] != 0.0:
                ratio = (-1) * U[k][k] / U[row][k]
                U[row] = [ui + uk / ratio for ui, uk in zip(U[row], U[k])]

    if permute_l:
        L = mat_mul(A, inverse(U))
    else:
        L = mat_mul(mat_mul(inverse(P), A), inverse(U))
        return (P, L, U)


def forward_sub(L, b):

    n = len(L)
    x = zeros(n)

    for i in range(n):
        tmp = b[i]
        for j in range(i):
            tmp -= L[i][j] * x[j]
        x[i] = tmp / L[i][i]
    
    return x


def back_sub(U, b):

    n = len(U)
    x = zeros(n)

    for i in range(n - 1, -1, -1):
        tmp = b[i]
        for j in range(i + 1, n):
            tmp -= U[i][j] * x[j]
        x[i] = tmp / U[i][i]

    return x


def solveLU(a, b):

    '''
        [A]{x} = {b} => [L][U]{x} = {b} => [L]{y} = {b}, {y} = [U]{x}
    '''

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


def solveLUP(A, b, tridiagonal=False):

    '''
        Solving system of equations with partial pivoting
        [A]{x} = {b} => [P][L][U]{x} = {b} => 
        => [P][L]{y} = {b}, [L]{y} = {b}(inverse([P])),
        where {y} = [U]{x} => {x} = {y}(inverse([U]))
    '''

    if tridiagonal:
        bm = [A[i][i] for i in range(n) for j in range(n) if i == j]
        am = [A[i+1][i] for i in range(n - 1) for j in range(n) if i == j]
        cm = [A[i][i+1] for i in range(n) for j in range(n - 1) if i == j]
        dm = deepcopy(b)
        x = zeros(n)

        for i in range(n - 1):
            w = am[i] / bm[i]
            bm[i + 1] -= w * cm[i]
            dm[i + 1] -= w * dm[i]
        x[n - 1] = dm[n - 1] / bm[n - 1]
        for i in range(n - 2, -1, -1):
            x[i] = (dm[i] - cm[i] * x[i + 1]) / bm[i]
        
        return x

    n = len(A)
    y = [0.0 for i in range(n)]
    x = [0.0 for i in range(n)]

    P, L, U = lup(A, permute_l=False)
    
    bP = mat_mul(b, inverse(P))

    # Forward substitution [L]{y} = {b}, {y} = [U]{x}
    for i in range(n):
        y[i] = (1 / L[i][i]) * (bP[i] - sum(L[i][j] * y[j] for j in range(i)))
    
    # Backward substitution [L]{x} = {y}
    for i in range(n-1, -1, -1):
        x[i] = (1 / U[i][i]) * (y[i] - sum(U[i][j] * x[j] for j in range(i+1, n)))
    
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


def jacobi(a: list, b: list, xinit: list = [0.0, 0.0, 0.0], es: float = 10 ** -6, iter: int = 100):
    '''
        Instead of solving [A]{x} = {b}, solve {x} = [C] + [M]{x}
            where [C] = inverse([D]){b}, [M] = inverse([D])[R],
            where [D] is a matrix [A] with all elements equal zero except of main diagonal elements 
            [R] is a matrix [A] with all elements of main diagonal equal zero
    '''

    if iter == 0:
        raise RecursionError('Maximum iterations exceeded. No solution found.')

    d = identity(len(a))
    r = deepcopy(a)
    for i in range(len(d)):
        d[i][i] = a[i][i]
        r[i][i] = 0.0
    
    d = inverse(d, True)
    c = mat_mul(d, b)
    m = mat_mul(d, r)

    x = mat_add(c, scalar_mul((-1), mat_mul(m, xinit)))

    if (er := (distance(x, xinit) / norm(x))) <= es:
        return x, er
    else:
        return jacobi(a, b, x, es, iter-1)


def gauss_seidel(a: list, b: list, xinit, es: float = 10 ** -6, iter: int = 100, w: float = 1.0):

    '''
        Instead of solving [A]{x}_(k+1) = {b}, solve [L]{x} = {b} - [U]{x}_(k)
            where [L] is lower diagonal matrix with elements of matrix [A]
            and [U] is upper diagonal matrix with elements of matrix [A] uii = 0
            w is a relaxation factor which is 0 < w < 2
    '''

    if not iter:
        raise RecursionError('Maximum iteration excceded. No solution found.')

    width = len(a)
    x = [0.0 for i in range(width)]

    for i in range(width):
        aii = 1 / a[i][i]
        s1 = sum([a[i][j] * xinit[j] for j in range(i + 1, width)])
        s2 = sum([a[i][j] * x[j] for j in range(0, i)])
        x[i] = aii * (b[i] - s1 - s2)
        x[i] = xinit[i] + w * (x[i] - xinit[i])

    if (er := (distance(x, xinit) / norm(x))) <= es:
        return x, er
    else:
        return gauss_seidel(a, b, x, es, iter-1)


# Solving Nonlinear systems of equations
def fixed_point(fx: list, x0: list, es: float = 10 ** -6, max_it: int = 100, current_it: int = 0):

    '''
        Solving nonlinear system of equations by using fixed point method
            fx is a list of funcions fx = [f1(x1, x2, ... xn), f2(x1, x2, ... xn), ..., fn(x1, x2, ... xn)]
            x0 is an initial guess. len(fx) must be equal len(x0)
            w is a relacsation factor which is 0 < w < 2
    '''

    x = [0.0 for i in range(len(fx))]
    tmp = [xi for xi in x0]
    if current_it >= max_it:
        raise RecursionError('Maximum iterations esceeded. No solutions found.')

    for i, fxi in enumerate(fx):
        x[i] = fxi(tmp)
        tmp[i] = x[i]

    er = abs(distance(x, x0) / norm(x))

    if er <= es:
        return x
    else:
        return fixed_point(fx, x, es, max_it, current_it + 1)


def newton_raphson(fx: list, jac: callable = None, x0: list = None, es = 10 ** -6, w: float | List[float] = None, iter: int = 0):

    '''
        Newton-Raphson method for solving nonlinear system of equations
            {f} = 
                f1(x1, x2, ..., xn) = 0
                f2(x1, x2, ..., xn) = 0
                ...
                fn(x1, x2, ..., xn) = 0

            {x(i+1)} = {x(i)} + {del_x}
            where {del_x} = -inverse([J]){f}
                  [J] is a Jacobian of {f} (nxn matrix)

            fx: List[callable]: list of equations
            jac: callable (default None): function that returns Jacobian matrix jac(*x0). If None Jacobian matrix wil
                    be computed numerically
            x0: list[float] (default None) is an initial guess, if None then will be set to 1.0 for all arguments
            w: float or list[float] (default None) is a relaxation factor 0.0 < w < 2.0 if None solve system without weight
    '''

    if iter > 800:
        raise RuntimeError('Maximum number of iterations exceeded. No solution found.')
    
    n = len(fx)
    nargs = len(signature(fx[0]).parameters)

    if n != nargs:
        raise RuntimeError('Undefined system of equations')

    if not x0: x0 = [1.0 for i in range(nargs)]

    if not jac:
        jac1 = [[0.0 for i in range(nargs)] for j in range(n)]
        for i in range(n):
            for j in range(nargs):
                dx = x[j] * 1e-6
                args1 = [x0[k] if k != j else x0[k] + dx for k in range(nargs)]
                args2 = [x0[k] if k != j else x0[k] - dx for k in range(nargs)]
                jac1[i][j] = (fx[i](*args1) - fx[i](*args2)) / (x * dx)
    else:
        jac1 = jac(*x0)

    func = [fxi(*x0) for fxi in fx]
    jac1 = inverse(jac1)
    delx = mat_mul(mat_mul((-1), jac1), func)
    if not w:
        x = mat_add(x0, delx)
    else:
        x = mat_add(x0, mat_mul(w, delx))
    er = norm(delx) / norm(x)

    if er <= es:
        return x
    else:
        return newton_raphson(fx, jac, x, es, w, iter + 1)


def newton_raphson2(f: callable, x0: float = None, args: list = [], tol: float = 1e-6, w: float = None ,it: int = 0):

    '''
        Solve equation of one unknown by using Newton-Raphson method
        ------------------------------------------------------------
        Parameters:     f: callable
                            Function of equation
                        x0: float (optional)
                            Initial guess. Update every iteration
                        args: array_like (optinal)
                            Extra argument to pass to the function.
                            Default is to pass no extra arguments
        Returns:        result: float
                            root of equation
    '''

    if it >= 800:
        raise RuntimeError('Maximum iteration exceeded. No solution found')
    
    if not x0 : x0 = 1.0
    
    der = (f(x0 + x0 * 1e-6, *args) - f(x0 - x0 * 1e-6, *args)) / (2 * x0 * 1e-6)
    xnew = x0 - f(x0, *args) / der if not w else x0 - w * f(x0, *args) / der
    es = abs(xnew - x0)
    if es <= tol:
        return xnew
    else:
        return newton_raphson2(f, xnew, args, tol, w, it  + 1)


def fsolve(func:callable, x0:list, args=(), dx:float=None, jac:callable=None, rtol:float=1e-6, max_iter:int=100, w:float=None) -> list[float]:
    
    '''
        Find roots of a function.
        Returns the roots of the (non-linear) equations defined by func(x) = 0 given a statring
        estimate
        ---------------------------------------------------------------------------------------
        Parameters:     
        -----------
                        func : callable f(x, *args)
                            A function that takes at least one argument, and returns a value 
                            of the same length
                        x0 : array_like, optional
                            The starting estimate for the roots of func(x) = 0.
                            If None will be set to 1.0 for all x
                        args : tuple, optional
                            Any extra arguments to func
                        dx : array_like, optional
                            Step size array. If None will be set to 1e-6 for each x
                        jac : callable, optional
                            Function to compute Jacobian matrix. If None will be estimated
                        rtol : float, optional
                            The calculation will terminate if the relative error between two
                            consecutive iterates is at most rtol. Default 1e-6
                        max_iter : integer, optional
                            Maximum number of iterations. Default 800
                        w : array_like, optional
                            Damping factor. If None estimated solution will not be damped
        Returns:        
        --------
                        x : array_like
                            The solution (or the result of last iteration for an unsuccessful call)
    '''

    
    n = len(x0)

    def jacobian_matrix(func, x, dx, args=()):
        jac = [[0.0 for j in range(n)] for i in range(n)]
        
        for i in range(n):
            for j in range(n):
                x1 = [x[k] if k != j else x[k] + dx[k] for k in range(n)]
                x2 = [x[k] if k != j else x[k] - dx[k] for k in range(n)]
                jac[i][j] = (func(x1, *args)[i] - func(x2, *args)[i]) / dx[j]
        return jac

    if not dx: dx = [1e-6 for i in range(n)]

    for i in range(max_iter):

        if not jac:
            jacobian = jacobian_matrix(func, x0, dx, args)
        else:
            jacobian = jac(x0, *args)
        
        fs = func(x0, *args)
        jacobian = inverse(jacobian)
        
        if not w:
            delx = mat_mul(mat_mul(-1, jacobian), fs)
        else:
            delx = mat_mul(mat_mul(mat_mul(-1, jacobian), fs), w)
        
        xn = mat_add(x0, delx)
        es = norm([(f0 - fn) ** 2 for f0, fn in zip(fs, func(xn, *args))]) / norm(fs)

        if abs(es) <= rtol:
            return xn
        else:
            x0 = xn

    print('Maximum iterration exceded. No solutions found')
    return xn
