import os
import sys


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(DIR)


from la.linalg import gauss_seidel as gs

def polyfit(points: list):

    n = len(points)
    x = [points[i][0] if i else 1.0 for i in range(n)]
    b = [points[i][1] for i in range(n)]
    a = [[x[j] ** i for i in range(n)] for j in range(n)]
    return gs(a, b, [0.0 for i in range(len(a))], iter=900)[0]