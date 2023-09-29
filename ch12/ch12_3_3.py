'''
    Solving ODE IVP by using Heun's method
'''

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__), os.pardir
))
sys.path.append(DIR)

from de.ode import euler, mid_point, heun, trapz,solve_ipv as sivp


def f1(t, x):
    return t * x ** 2 + 2 * x


def f1exact(t, x0):
    return (-20 * np.exp(2 * t)) / (9 - 5 * np.exp(2 * t) + 10 * np.exp(2 * t) * t)


def vdp(t, y, mu):
    x = y[0]
    v = y[1]
    return [v, mu * (1 - x ** 2) * v - x]


def f2(t, x, omega):
    y1 = x[0]
    y2 = x[1]
    return [y2, (-omega ** 2) * y1]


def f2exact(t, y0, omega):
    return y0 * np.cos(omega * t)


def f3(t, y):
    return y ** 2 - y


def f3ex(t):
    return 1 / (1 + np.exp(t))


if __name__ == "__main__":

    y0 = -5.0
    tspan = [0, 5]
    h = 0.1
    teval = np.linspace(tspan[0], tspan[1], 50)
    teu, yeu = euler(f1, tspan, y0, h)
    theun, yheun = heun(f1, tspan, y0, h)
    tmid, ymid = mid_point(f1, tspan, y0, h)
    trk, yrk = sivp(f1, tspan=tspan, ic=y0, teval=teval)

    tex = np.linspace(0, 5, 200)
    xex = f1exact(tex, y0)

    fig, ax = plt.subplots()
    ax.plot(tex, xex, marker='', color='k', label='exact')
    ax.plot(teu, yeu, marker='', linestyle='-.', label='Euler')
    ax.scatter(theun, yheun, marker='o', label='Heun')
    ax.scatter(tmid, ymid, marker='*', label='Mid Point')
    ax.scatter(trk, yrk, marker='^', color='r', label='RK45')

    
    plt.legend()
    plt.show()

    t = np.linspace(0, 6, 40)
    tex = np.linspace(0, 6, 600)
    omega = 4
    yex = f2exact(tex, 1, 4)
    teu, yeu = euler(f2, t_span=[t[0], t[-1]], y0=[0, 1], h=0.15, args=[omega])
    th, yh = heun(f2, t_span=[t[0], t[-1]], y0=[0, 1], h=0.15, args=[omega])
    tm, ym = mid_point(f2, t_span=[t[0], t[-1]], y0=[0, 1], h=0.15, args=[omega])
    ysol = solve_ivp(f2, t_span=[t[0], t[-1]], y0=[0, 1], args=[omega], t_eval=t)
    ts, ys = sivp(f2, [t[0], t[-1]], [0, 1], args=[omega])
    tt, yt = trapz(f2, [t[0], t[-1]], [0, 1], 0.15, [omega])
    yeu1 = np.array(yeu)[:, 1]
    yh1 = np.array(yh)[:, 1]
    ym1 = np.array(ym)[:, 1]
    yt1 = np.array(yt)[:, 1]

    fig, ax = plt.subplots()
    ax.plot(tex, yex, marker='', label='Exact', color='k')
    ax.plot(teu, yeu1, marker='o', color='r', label='Euler')
    ax.scatter(th, yh1, marker='<', color='b', label='Heun')
    ax.scatter(tm, ym1, marker='*', color='g', label='MidPoint')
    ax.scatter(ysol.t, ysol.y[1], marker='s', color='k', label='Scipy')
    ax.scatter(ts, ys[1], marker='o', color='g', label='RK45')
    ax.plot(tt, yt1, marker='*', linestyle='-.', color='k', label='Trapz')

    ax.set_ylim(-1.5, 1.5)


    plt.legend()
    plt.show()

    mu = 0
    t = np.linspace(0, 50, 500)
    tsol, ysol = heun(vdp, t_span=[t[0], t[-1]], y0=[1, 0], h=0.1, args=[mu])
    sol = solve_ivp(vdp, t_span=[t[0], t[-1]], y0=[1, 0], t_eval=t, args=[mu])
    tt, yt = trapz(vdp, t_span=[t[0], t[-1]], y0=[1, 0], h=0.1, args=[mu])

    ysol1, ysol2 = np.array(ysol)[:, 0], np.array(ysol)[:, 1]
    yt1, yt2 = np.array(yt)[:, 0], np.array(yt)[:, 1]

    fig, ax = plt.subplots()
    ax.scatter(tsol, ysol2, marker='o', color='r')
    ax.scatter(tt, yt2, marker='^', color='y')
    ax.plot(sol.t, sol.y[1])

    plt.show()

    y0 = 0.5
    t = np.linspace(0, 1, 10)
    tex = np.linspace(0, 1, 100)
    yex = f3ex(tex)
    tsol, ysol = trapz(f3, [t[0], t[-1]], 0.5, 0.1)

    fig, ax = plt.subplots()
    ax.plot(tex, yex, marker='', color='k')
    ax.scatter(tsol, ysol, marker='o', color='b')

    plt.show()
