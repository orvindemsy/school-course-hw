'''
This code will include three ways to apply finite difference approximation, which are:
- Gauss-Seidel (GS)
- Gauss Elimination (GE)
- Conjugate Gradient (CG)
'''

from numpy import *
from matplotlib.pyplot import *
import time

# definition of computational domain
w = 10.0
h = 15.0

# Conjugate Gradient
def conjgrad(a, b, eps):
    n = b.size
    imax = 10*n
    f = zeros(n, float)
    r0 = b - dot(a,f)
    p = r0
    for i in range(0, imax):
        alpha = dot(r0, r0)/dot(p, dot(a, p))
        f = f + alpha*p
        r1 = r0 - alpha*dot(a, p)
        res = linalg.norm(r1)
        beta = dot(r1, r1)/dot(r0, r0)
        p = r1 + beta*p
        r0 = r1
        # print (i, res)
        if res < eps: return f
    # print ("not converged!")
    return f

# Gauss Elimination
def gausselim(a, b):
    n = b.size
    f = zeros(n, float)
    for k in range(0,n-1):
        for i in range(k+1,n):
            em = a[i,k]/a[k,k]
            a[i,k] = em
            b[i] -= em*b[k]
            for j in range(k+1,n):
                a[i,j] -= em*a[k,j]

    f[n-1] = b[n-1]/a[n-1,n-1]

    for i in range(n-2,-1,-1):
        f[i] = b[i]
        for j in range(n-1,i,-1):
            f[i] -= a[i,j]*f[j]
        f[i] = f[i]/a[i,i]

    return f

# Gauss Seidel
def gauss_seidel(T, nx, ny, beta, eps):
    itrmax = 1000
    bsq = beta*beta
    dnm = 0.5/(1.0+bsq)
    for itr in range(0, itrmax):
        err = 0.0
        for j in range(1, ny):
            for i in range(1, nx):
                new = (T[i+1,j]+bsq*T[i,j+1]+T[i-1,j]+bsq*T[i,j-1])*dnm
                err = max(abs(T[i,j]-new),err)
                T[i,j] = new
        if (err < eps): return itr
    return itrmax


def time_elapsed(t0):
    t1 = time.time()
    return t1-t0


def laplace_CG(nx, ny):
    start = time.time()
    nmax = (nx - 1) * (ny - 1)
    dx = w / nx
    dy = h / ny
    beta = dx / dy
    bsq = beta * beta
    dnm = 0.5 / (1.0 + bsq)

    # initialization
    x = linspace(0, w, nx + 1)
    y = linspace(0, h, ny + 1)
    b = zeros(nmax, float)
    a = zeros((nmax, nmax), float)
    T = zeros((nx + 1, ny + 1), float)

    # Dirichlet boundary condition for upper boundary (j = ny)
    for i in range(0, nx + 1):
        T[i, ny] = 100.0 * sin(pi * x[i] / w)

    n = 0
    for j in range(1, ny):
        for i in range(1, nx):
            if (j == 1):
                b[n] += bsq * dnm * T[i, j - 1]  # lower boundary
            else:
                a[n - nx + 1, n] = -bsq * dnm
            if (i == 1):
                b[n] += dnm * T[i - 1, j]  # left boundary
            else:
                a[n - 1, n] = -dnm
            if (i == nx - 1):
                b[n] += dnm * T[i + 1, j]  # right boundary
            else:
                a[n + 1, n] = -dnm
            if (j == ny - 1):
                b[n] += bsq * dnm * T[i, j + 1]  # upper boundary
            else:
                a[n + nx - 1, n] = -bsq * dnm
            a[n, n] = 1.0
            n += 1

    eps = 1e-8
    f = conjgrad(a, b, eps)
    n = 0
    for j in range(1, ny):
        for i in range(1, nx):
            T[i, j] = f[n]
            n += 1

    # print("Time for CG : ", time_elapsed(start))
    #
    # X, Y = meshgrid(x, y, indexing='ij')
    # contour(X, Y, T)
    # xlabel('X')
    # ylabel('Y')
    # cbar = colorbar()
    # cbar.set_label('T')
    # show()

    return time_elapsed(start)


def laplace_GE(nx, ny):
    start = time.time()
    nmax = (nx - 1) * (ny - 1)
    dx = w / nx
    dy = h / ny
    beta = dx / dy
    bsq = beta * beta
    dnm = 0.5 / (1.0 + bsq)

    x = linspace(0, w, nx + 1)
    y = linspace(0, h, ny + 1)
    b = zeros(nmax, float)
    a = zeros((nmax, nmax), float)
    T = zeros((nx + 1, ny + 1), float)

    # Dirichlet boundary condition for upper boundary (j = ny)
    for i in range(0, nx + 1):
        T[i, ny] = 100.0 * sin(pi * x[i] / w)

    n = 0
    for j in range(1, ny):
        for i in range(1, nx):
            if (j == 1):
                b[n] += bsq * dnm * T[i, j - 1]  # lower boundary
            else:
                a[n - nx + 1, n] = -bsq * dnm
            if (i == 1):
                b[n] += dnm * T[i - 1, j]  # left boundary
            else:
                a[n - 1, n] = -dnm
            if (i == nx - 1):
                b[n] += dnm * T[i + 1, j]  # right boundary
            else:
                a[n + 1, n] = -dnm
            if (j == ny - 1):
                b[n] += bsq * dnm * T[i, j + 1]  # upper boundary
            else:
                a[n + nx - 1, n] = -bsq * dnm
            a[n, n] = 1.0
            n += 1

    f = gausselim(a, b)
    n = 0
    for j in range(1, ny):
        for i in range(1, nx):
            T[i, j] = f[n]
            n += 1

    # print(T)
    # print("Time for GE : ", time_elapsed(start))

    # X, Y = meshgrid(x, y, indexing='ij')
    # contour(X, Y, T)
    # xlabel('X')
    # ylabel('Y')
    # cbar = colorbar()
    # cbar.set_label('T')
    # show()

    return time_elapsed(start)


def laplace_GS(nx, ny):
    start = time.time()
    dx = w / nx
    dy = h / ny
    beta = dx / dy

    # initialization
    x = linspace(0, w, nx + 1)
    y = linspace(0, h, ny + 1)
    T = zeros((nx + 1, ny + 1), float)

    # Dirichlet boundary condition for upper boundary (j = ny)
    for i in range(0, nx + 1):
        T[i, ny] = 100.0 * sin(pi * x[i] / w)

    eps = 1e-8
    itr = gauss_seidel(T, nx, ny, beta, eps)

    # obtained solution in standard output
    # time_e = time_elapsed(start)
    # print("Time for GS : ", time_elapsed(start))

    # # plot setting
    # X, Y = meshgrid(x, y, indexing='ij')
    # contour(X, Y, T)
    # xlabel('X')
    # ylabel('Y')
    # cbar = colorbar()
    # cbar.set_label('T')
    # show()

    return time_elapsed(start)