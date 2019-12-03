from numpy import *
from matplotlib.pyplot import *
import time

# GE stands for grad elimination
# Direct method
# b is the right hand side so should be defined
# bsq is b square
def gausselim(a,b):
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

# For direct method we have to define the matrix this is the definition part of the matrix
# Page 539 in textbook
# n, n is the diagonal element
# In the textbook it's -4 but it's just normalization technique
# b is defined at right hand side so its important
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

f = gausselim(a,b)
n = 0
for j in range(1, ny):
    for i in range(1, nx):
        T[i,j] = f[n]
        n += 1

# obtained solution in standard output
time_e = time_elapsed(start)
print ("solution =", T)
print ("elapsed time =",time.time()-start)
# print ("elapsed time =", time_e)

# print(T)

# plot setting
X, Y = meshgrid(x,y,indexing='ij')
contour(X,Y,T)
xlabel('X')
ylabel('Y')
cbar=colorbar()
cbar.set_label('T')
show()

# Needs larger operation for larger data
# because it's proportional to cubic operation