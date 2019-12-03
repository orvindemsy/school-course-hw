from numpy import *
from matplotlib.pyplot import *
import time

# The definition is the same as GE or GS
# CG method is very efficient

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
        print (i, res)
        if res < eps: return f
    print ("not converged!")
    return f

# Added time elapsed function
def time_elapsed(t0):
    t1 = time.time()
    return t1-t0

start = time.time()

# definition of computational domain
w = 10.0
h = 15.0
nx = 10
ny = 10
nmax = (nx-1)*(ny-1)
dx = w/nx
dy = h/ny
beta = dx/dy
bsq = beta*beta
dnm = 0.5/(1.0+bsq)

# initialization
x = linspace(0, w, nx+1)
y = linspace(0, h, ny+1)
b = zeros(nmax, float)
a = zeros((nmax,nmax), float)
T = zeros((nx+1,ny+1), float)

# Dirichlet boundary condition for upper boundary (j = ny)
for i in range(0, nx+1):
    T[i,ny] = 100.0*sin(pi*x[i]/w)

n = 0
for j in range(1, ny):
    for i in range(1, nx):
        if ( j==1 ):
            b[n] += bsq*dnm*T[i,j-1]      # lower boundary
        else:
            a[n-nx+1,n] = -bsq*dnm
        if ( i==1):
            b[n] += dnm*T[i-1,j]          # left boundary
        else:
            a[n-1,n] = -dnm
        if ( i==nx-1 ):
            b[n] += dnm*T[i+1,j]          # right boundary
        else:
            a[n+1,n] = -dnm
        if ( j==ny-1 ):
            b[n] += bsq*dnm*T[i,j+1]      # upper boundary
        else:
            a[n+nx-1,n] = -bsq*dnm
        a[n,n] = 1.0
        n += 1

eps = 1e-8
f = conjgrad(a,b,eps)
n = 0
for j in range(1, ny):
    for i in range(1, nx):
        T[i,j] = f[n]
        n += 1

# obtained solution in standard output
print("solution =", T)

time_e = time_elapsed(start)
# print ("elapsed time =",time.time()-start)
print("elapsed time =", time_e)

# plot setting
X, Y = meshgrid(x, y, indexing='ij')
contour(X,Y,T)
xlabel('X')
ylabel('Y')
cbar=colorbar()
cbar.set_label('T')
show()
