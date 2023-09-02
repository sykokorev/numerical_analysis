def polyint(x, *coeff):
    return sum([(1 / (i + 1)) * c * x ** (i + 1) for i, c in enumerate(coeff)])
