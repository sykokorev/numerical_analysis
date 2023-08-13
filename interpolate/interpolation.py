import os
import sys


DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    os.pardir
))
sys.path.append(DIR)


from la.linalg import *


def polyfit(points: list):

    '''
        Function returns coefficients of the polynomial
        which fits the points
        The points are list of points [[x1, y1], ... [xn, yn]]
    '''

    n = len(points)
    x_p = [points[i][0] if i else 1.0 for i in range(n)]
    x = []
    
    for i, p in enumerate(x_p):
        tmp = [1 if j == 0 else p ** j for j in range(n)]
        x.append(tmp)

    b = [points[i][1] for i in range(n)]
    x = inverse(x)
    a = mat_mul(b, x)
    return a


def divided_diff(x, y, b=[]):
    '''
        Function to calculate the divided
        difference
    '''
    n = len(x)

    c = [[0.0 for i in range(n)] for j in range(n)]
    for i, yi in enumerate(y):
        c[i][0] = yi

    for j in range(1, n):
        for i in range(n - j):
            c[i][j] = (c[i+1][j-1] - c[i][j-1]) / (x[i+j] - x[i])
    
    return c[0][:]


def Newtoninterp(x, y, xp):

    '''
        Evaluate the newton polynomial at xp
    '''
    n = len(x) - 1
    c = divided_diff(x, y)
    yp = c[n]
    for k in range(1, n+1):
        yp = c[n-k] + (xp - x[n-k]) * yp
    return yp


def Lagrangeinterp(x, y, xp):

    n = len(x)
    yp = 0.0
    for i in range(n):
        L = 1.0
        for j in range(n):
            if j != i:
                L *= (xp - x[j]) / (x[i] - x[j])
        yp += L * y[i]
    return yp
