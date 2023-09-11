import os
import json
import numpy as np

from integration.integrate import *
from la.linalg import *



FILE = os.path.abspath(os.path.join(
    os.path.dirname(__file__), 
    'integration',
    'guassian_quadrature.json'
))


if __name__ == "__main__":

    points = []
    gauss_quad_json = []

    for n in range(1, 40):

        num_points = n
        args = legendre_poly_coef(num_points)
        intervals = np.linspace(-1, 1, 300)

        xs = []
        ws = []

        s = sum([a[0] for a in args])
        args_norm = [(a[0] / s, a[1]) for a in args]

        for i, it in enumerate(intervals):
            try:
                xi = newton_raphson2(legendre_poly, it, args)
                der = legendre_der(xi, *args_norm)
                if i:
                    if all([abs(xi - xs[j]) > 1e-3 for j in range(len(xs))]): 
                        xs.append(xi)
                        ws.append(2 / ((1 - xi ** 2) * (der ** 2)))
                else:
                    xs.append(xi)
                    ws.append(2 / ((1 - xi ** 2) * (der ** 2)))
            except RuntimeError:
                print(f'{n}, {i}. Root no found')
        points.append((round(xi, 8), round(wi, 8)) for xi, wi in zip(xs, ws))
        gauss_quad_json.append(list(points[n - 1]))

    with open(FILE, 'w') as fp:
        json.dump(gauss_quad_json, fp)
