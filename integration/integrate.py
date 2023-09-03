def polyint(x, *coeff):
    return sum([(1 / (i + 1)) * c * x ** (i + 1) for i, c in enumerate(coeff)])


def rectangle(f, x, args = [], method = 'mid'):

    I = 0
    if method == 'mid':
        I = sum([f((x[i] + x[i - 1]) / 2) * (x[i] - x[i - 1]) for i in range(1, len(x))])
    elif method == 'right':
        I = sum([f(x[i - 1]) * (x[i] - x[i - 1]) for i in range(1, len(x))])
    elif method == 'left':
        I = sum([f(x[i]) * (x[i] - x[i - 1]) for i in range(1, len(x))])

    return I


def trapezoid(f, x, args = []):

    n = len(x)
    I = 0.0
    for i in range(1, n):
        h = x[i] - x[i - 1]
        I += h * (f(x[i - 1], *args) + f(x[i], *args)) / 2
    
    return I


def simpson(f, x, args = [], method=1):

    '''
        method = 1 is Simpson's 1/3 rule
        method = 2 is Simpson's 3/8 rule
    '''

    n = len(x)
    I = 0.0
    if method == 1:
        for i in range(1, n):
            h = (x[i] - x[i - 1]) / 2
            xm = x[i - 1] + h
            I += (f(x[i - 1], *args) + 4 * f(xm, *args) + f(x[i], *args)) * (h / 3)
    if method == 2:
        for i in range(1, n):
            h = (x[i] - x[i - 1]) / 3
            xli = x[i - 1] + h
            xri = x[i - 1] + 2 * h
            I += (f(x[i - 1], *args) + 3 * f(xli, *args) + 3 * f(xri, *args) + f(x[i])) * (3 * h / 8)

    return I
