def derivative(f, x: float, n: int = 1, args: list = [], dx: float = 1e-4   ):

    if n == 1:
        return ((-1) * f(x + 2 * dx, *args) + 8 * f(x + dx, *args) - 8 * f(x - dx, *args) + f(x - 2 * dx, *args)) / (12 * dx ** n)
    elif n == 2:
        return ((-1) * f(x + 2 * dx, *args) + 16 * f(x + dx, *args) -30 * f(x, *args) + 16 * f(x - dx, *args) - f(x - 2 * dx, *args)) / (12 * dx ** n)
    elif n == 3:
        return ((-1) * f(x + 3 * dx, *args) + 8 * f(x + 2 * dx, *args) - 13 * f(x + dx, *args) + 13 * f(x - dx, *args) - 8 * f(x - 2 * dx, *args) + f(x - 3 * dx, *args)) / (8 * dx ** n)
    elif n == 4:
        return ((-1) * f(x + 3 * dx, *args) + 12 * f(x + 2 * dx, *args) - 39 * f(x + dx, *args) + 56 * f(x, *args) - 39 * f(x - dx, *args) + 12 * f(x - 2* dx, *args) - f(x - 3 * dx, *args)) / (6 * dx ** n)
    else:
        return None


def differnce(xdata: list, ydata: list, n: int = 1):

    k = len(xdata) - 1
    y = []
    x = []
    h = xdata[1] - xdata[0]

    if n == 1:
    
        for i in range(2, k - 2):
            x.append(xdata[i])
            y.append(((-1) * ydata[i + 2] + 8 * ydata[i + 1] - 8 * ydata[i - 1] + ydata[i - 2]) / (12 * h))
    
    elif n == 2:

        for i in range(2, k - 2):
            x.append(xdata[i])
            y.append(
                ((-1) * ydata[i + 2] + 16 * ydata[i + 1] -30 * ydata[i] + 16 * ydata[i - 1] - ydata[i - 2]) / (12 * h ** n)
            )
    
    elif n == 3:

        for i in range(3, k - 3):
            x.append(xdata[i])
            y.append(
                ((-1) * ydata[i + 3] + 8 * ydata[i + 2] - 13 * ydata[i + 1] + 13 * ydata[i - 1] - 8 * ydata[i - 2] + ydata[i - 3]) / (8 * h ** n)
            )
    
    elif n == 4:

        for i in range(3, k - 3):
            x.append(xdata[i])
            y.append(
                ((-1) * ydata[i + 3] + 12 * ydata[i + 2] - 39 * ydata[i + 1] + 56 * ydata[i] - 39 * ydata[i - 1] + 12 * ydata[i - 2] - ydata[i - 3]) / (6 * h ** n)
            )

    return (x, y)
