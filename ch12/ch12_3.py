'''
    Solving ODE
    -----------
    Solving IVPs with using (explicit) Euler method
'''

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import math


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__), os.pardir
))
sys.path.append(DIR)

from ode.euler import *


'''
    Population grown rate problem
'''
def population(t: float, x: float, *args) -> float:
    return args[0] * x


def exact_pop(x0, x, k):
    return x0 * np.exp(k * x)


def concentration(t, c, *args):
    return (100 + 50 * np.cos(2 * np.pi * t / 365)) * (5 * np.exp((-1) * 2 * t / 1000) - c) / args[0]


def f2(t, x):
    return -0.5 * x


def f2_exact(x0, x):
    return x0 * np.exp(-0.5 * x)


if __name__ == "__main__":

    '''
        x represents population in million
        t represents time in years
            Initial condition x0 = 35 million, t0 = 0
            step h = ti+1 - ti = 1 year
    '''
    ts, ys = euler(population, t_span=(0, 50), y0=35, h=1.2, args=[0.015])
    exact = exact_pop(35, 50, 0.015)
    print(f'Approximated solution = {ys[-1]}')
    print(f'Exact solution = {exact}')
    e = abs(exact - ys[-1])
    er = e / exact
    print(f'E = {e}')
    print(f'Er = {er}')
    fig, ax = plt.subplots()
    ax.plot(ts, ys, marker='o', color='r')
    xex = np.linspace(0, 50, 50)
    yex = exact_pop(35, xex, 0.015)
    ax.plot(xex, yex, marker='', color='k')

    plt.show()

    '''
        Given a lake V = 10000 m^3
        Initial concentration of pollutant c(t0) = 5 (ppb). The soil is no longer polluted.
        Flow rate in and out of the lake q = 100 + 50 cos(2 pi t / 365) [m^3 / day]
        Decreasing concentration of pollutant in the surrounding soil cin = 5e^(-2 t / 1000)
        Total amount the pollutant in the lake each day = [100 + 50 cos(2 pi t / 365)] (cin - c)
        Concntration = [100 + 50 cos(2 pi t / 365)] (cin - c) / 10000
        Find concentration after 2 years
    '''
    volume = 10000
    c0 = 5
    t0 = 0
    tf = 2 * 365
    ts, ys = euler(concentration, t_span=(t0, tf), y0=c0, h=1, args=[volume])

    print(f'Solution of the task of the lake pollution is {ys[-1]}')

    fig, ax = plt.subplots()
    ax.plot(ts, ys, marker='', linestyle='--', color='b')
    plt.show()

    print(f'Exact solution of the function f2 equals {f2_exact(1, 20)}')
    step = 1
    ts, ys = euler(f2, (0, 20), 1, step)
    print(f'Approximated solution of the function f2 with step {step} equals {ys[-1]}')
    step2 = 4.2
    ts2, ys2 = euler(f2, (0, 20), 1, step2)
    print(f'Approximated solution of the function f2 with step {step2} equals {ys2[-1]}')
    xex = np.linspace(0, 20, 1000)
    fex = f2_exact(1, xex)

    fog, ax = plt.subplots()
    ax.plot(xex, fex, marker='', linestyle='-', color='k', label='Exact')
    ax.plot(ts, ys, marker='', linestyle='--', color='r', label=f'Step-size {step}')
    ax.plot(ts2, ys2, marker='', linestyle='-.', color='g', label=f'Step-size {step2}')

    plt.legend()
    plt.show()
