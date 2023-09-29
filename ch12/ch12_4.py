import os
import sys
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__), os.pardir
))
sys.path.append(DIR)


from de.ode import solve_ipv as sivp, euler


def spring_mass_system(t, x, m, k, c):

    '''
        m x'' + c x' + k x = 0
            x' = y
            y' = - (c y + k x) / m
    '''

    x1 = x[0]   # dependet variable x(t)
    x2 = x[1]   # x'(t)
    return [
        x2,
        (-1) * (c * x2 + k * x1) / m
    ]


if __name__ == "__main__":

    x0 = [-10, 0]
    tspan = [0, 20]
    h = 0.1
    t = np.linspace(0, 20, 200)

    args = [1.0, 1.0, 0.15]

    tsol, ysol = sivp(spring_mass_system, tspan, x0, teval=t, args=args)
    teu, yeu = euler(spring_mass_system, tspan, x0, h, args=args)
    yeu1, yeu2 = np.array(yeu)[:, 0], np.array(yeu)[:, 1]
    sol = sc.integrate.solve_ivp(spring_mass_system, tspan, x0, t_eval=t, args=args)


    fig, axs = plt.subplots(nrows=1, ncols=2)

    axs[0].plot(tsol, ysol[0], label='RK45')
    axs[0].plot(teu, yeu1, label='Explicit Euler', marker='', linestyle='--', color='r')
    axs[0].plot(sol.t, sol.y[0], label='SciPy RK45', marker='', linestyle='--', color='b')
    axs[1].plot(tsol, ysol[1], label='RK45')
    axs[1].plot(teu, yeu2, label='Explicit Euler', marker='', linestyle='--', color='r')
    axs[1].plot(sol.t, sol.y[1], label='SciPy RK45', marker='', linestyle='--', color='b')

    axs[0].set_title('Position')
    axs[0].set_xlabel('Time, [s]')
    axs[0].set_ylabel('Position, [m]')

    axs[1].set_title('Velocity')
    axs[1].set_xlabel('Time, [s]')
    axs[1].set_ylabel('Velocity, [m/s]')

    for ax in axs:
        ax.grid(True)

    plt.legend()
    plt.show()
