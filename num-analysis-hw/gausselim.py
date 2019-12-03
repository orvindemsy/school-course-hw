'''
Source: http://www.rhd.mech.tohoku.ac.jp/~ohnishi/free/gausselim.py
Example of a Gauss Elimination method
'''
from numpy import *


def gauss(a, b):
    n = b.size
    x = zeros(n, float)
    for k in range(0, n-1):
        for i in range(k+1, n):
            em = a[i, k]/a[k, k]
            a[i, k] = em
            b[i] -= em*b[k]
            for j in range(k+1, n):
                a[i, j] -= em*a[k, j]

    x[n-1] = b[n-1]/a[n-1, n-1]

    for i in range(n-2, -1, -1):
        x[i] = b[i]
        for j in range(n-1, i, -1):
            x[i] -= a[i, j]*x[j]
        x[i] = x[i]/a[i, i]

    return x


a = array([[80.0, -20.0, -20.0], [-20.0, 40.0, -20.0], [-20.0, -20.0, 130.0]])
b = array([20.0, 20.0, 20.0])

print(a, ",", b)
sol = gauss(a, b)
print(sol)
