def deriv1D(xi: list, yi: list):
    '''
        f: callable: Function to take partial derivative. 
            First and second arguments 
        xi and xj: list: grid points in the i- and j-directions
        args: list: Function arguments f(x, )
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

def deriv2D(xi: list, xj: list):
    
    n = len(xi)

    d = [0.0 for i in range(n)]
    for i in range(n):
        pass


def deriv3D(f, xi: list, xj: list, xk: list, args: list, h: float = 1e-3):
    pass
