def deriv1D(xi: list, yi: list):
    '''
        Compute difference 1D non-uniform grid
        xi and xj: list: grid points in the i- and j-directions
    '''
    
    n = len(xi)

    d = [0.0 for i in range(n)]
    for i in range(n):

        if i == 0:
            dx = xi[i + 1] - xi[i]
            d[i] = (yi[i + 1] - yi[i]) / dx
        elif i == n - 1:
            dx = xi[i] - xi[i - 1]
            d[i] = (yi[i] - yi[i - 1]) / dx
        else:
            dx = (xi[i + 1] - xi[i]) + (xi[i] - xi[i - 1])
            d[i] = (yi[i + 1] - yi[i - 1]) / dx

    return d

def deriv2D(xi: list, xj: list, xk):
    
    '''
        Compute difference 2D non-uniform grid
    '''

    n = len(xi)

    d = [0.0 for i in range(n)]
    for i in range(n):
        if i == 0:
            di = xi[i + 1] - xi[i]
            dj = xj[i + 1] - xj[i]
            d[i] = (1 / (di)) * ((xk[i + 1] - xk[i]) / dj)
        elif x == n - 1:
            di = xi[i] - xi[i - 1]
            dj = xj[i] - xj[i - 1]
            d[i] = (1 / (di)) * ((xk[i] - xk[i - 1]) / dj)
        else:
            di = (xi[i + 1] - xi[i]) + (xi[i] - xi[i - 1])
            dj = (xj[i + 1] - xj[i]) + (xj[i] - xj[i - 1])
            d[i] = (1 / di) * ((xk[i + 1] - xk[i - 1]) / dj)
        
    return d
