import math


def tiso(t, pi, t1, *args):

    const = args[0] * math.log(t1) + args[1] * t1 + (args[2] * t1 ** 2) / 2 + (args[3] * t1 ** 3) / 3 + (args[4] * t1 ** 4) / 4 + math.log(pi)
    return args[0] * math.log(t) + args[1] * t + (args[2] * t ** 2) / 2 + (args[3] * t ** 3) / 3 + (args[4] * t ** 4) / 4 - const


def find_t(t, pi, t1, args, es=1e-6, it = 0):

    if it >= 800:
        raise RuntimeError('Maximum number of iteration exceeded. No solution found.')

    dt = t * 1e-6
    df = (tiso(dt + t, pi, t1, *args) - tiso(t - dt, pi, t1, *args)) / (2 * dt)
    tnew = t - tiso(t, pi, t1, *args) / df

    if abs(tnew - t) <= es:
        return tnew
    else:
        return find_t(tnew, pi, t1, args, es, it + 1)


def const(t, *args):
    return args[0] * t + (args[1] * t ** 2) / 2 + (args[2] * t ** 3) / 3 + (args[3] * t ** 4) / 4 + (args[4] * t ** 5) / 5

def eff(t1, t2iso, t2, *args):

    c1 = const(t1, *args)
    c2 = const(t2, *args)
    c2iso = const(t2iso, *args)
    return (c2iso - c1) / (c2 - c1)


if __name__ == "__main__":

    a = [3.5683962, -0.000678729429, 0.00000155371476, -3.2993706e-12, -4.66395387e-13]
    # Axial
    pr = [1.247819, 1.2880493, 1.3296493]
    
    t2 = [366.22, 370.19, 374.38]

    t1 = 288.15
    t2iso = []
    effs = []

    for pi, t2i in zip(pr, t2):
        t = find_t(t2i, pi, t1, a, es=0.1)
        t2iso.append(t)
        effs.append(eff(t1, t, t2i, *a))

    print('Axial')
    print(str([round(effi, 4) for effi in effs])[1:-1])
    print(str([round(ti, 3) for ti in t2iso])[1:-1])

    t2iso = []
    effs = []

    pr = [3.745, 4.171, 4.749, 5.420, 5.088, 5.150, 5.184, 5.226, 4.940]
    t2 = [469.430, 480.024, 488.734, 498.077, 488.947, 489.907, 490.787, 491.887, 486.138]

    for pi, t2i in zip(pr, t2):
        t = find_t(t2i, pi, t1, a, es=0.1)
        t2iso.append(t)
        effs.append(eff(t1, t, t2i, *a))

    print('Centrigugal')
    print(str([round(effi, 4) for effi in effs])[1:-1])
    print(str([round(ti, 3) for ti in t2iso])[1:-1])
