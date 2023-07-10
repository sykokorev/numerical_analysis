import matplotlib.pyplot as plt
import numpy as np


def fact(n: int) -> int:

    if n <= 1:
        return 1
    else:
        return n * fact(n - 1)


def cos(x: float, n: int) -> float:

    if n < 1:
        return 1

    return ((-1) ** n) * (x ** (2 * n)) / fact(2 * n) + cos(x, n - 1)


def sin(x: float, n: int) -> float:

    if n < 1:
        return x

    return ((-1) ** n) * (x ** (2 * n + 1)) / fact(2 * n + 1) + sin(x, n - 1)


def exp(x: float, n: int) -> float:

    if n < 1:
        return 1
    
    return (x ** n) / fact(n) + exp(x, n - 1)



if __name__ == "__main__":

    x = 10
    n = np.linspace(0, 80, 81)
    c = cos(x, n[0])
    er = [0]
    cs = [c]
    es = 10 ** -6

    for i, o in enumerate(n[1:]):
        ci = cos(x, o)
        cs.append(ci)
        er.append((ci - c) / ci)
        if er[i + 1] < es: break
        c = ci
    
    fig, axs = plt.subplots(2, 1)
    axs[0].plot(n[:i+2], cs)
    axs[1].plot(n[:i+2], er)
    axs[0].grid(True)
    axs[0].set_title('Cos')
    axs[1].grid(True)
    axs[1].set_title('Error')

    plt.show()
