def diff(f, x: float, n: int = 1, args: list = [], dx: float = 1e-1):

    dx *= x

    if n == 1:
        return ((-1) * f(x + 2 * dx, *args) + 8 * f(x + dx, *args) - 8 * f(x - dx, *args) + f(x - 2 * dx, *args)) / (12 * dx ** n)
    elif n == 2:
        return ((-1) * f(x + 2 * dx, *args) + 16 * f(x + dx, *args) -30 * f(x, *args) + 16 * f(x - dx, *args) - f(x - 2 * dx, *args)) / (12 * dx ** n)
    elif n == 3:
        return ((-1) * f(x + 3 * dx, *args) + 8 * f(x + 2 * dx, *args) - 13 * f(x + dx, *args) + 13 * f(x - dx, *args) - 8 * f(x - 2 * dx, *args) + f(x - 3 * dx, *args)) / (8 * dx ** n)
    elif n == 4:
        return ((-1) * f(x + 3 * dx, *args) + 12 * f(x + 2 * dx, *args) - 39 * f(x + dx, *args) + 56 * f(x, *args) - 39 * f(x - dx, *args) + 12 * f(x - 2* dx, *args) - f(x - 3 * dx, *args)) / (6 * dx ** n)
