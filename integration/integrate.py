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
