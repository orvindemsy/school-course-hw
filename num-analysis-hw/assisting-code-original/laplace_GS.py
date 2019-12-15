from numpy import *
from matplotlib.pyplot import *
import time

def gauss_seidel(T,nx,ny,beta,eps):
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

start = time.time()

# definition of computational domain
w = 10.0
h = 15.0
nx = 10
ny = 10
dx = w/nx
dy = h/ny
beta = dx/dy

# initialization
x = linspace(0, w, nx+1)
y = linspace(0, h, ny+1)
T = zeros((nx+1,ny+1), float)

# Dirichlet boundary condition for upper boundary (j = ny)
for i in range(0, nx+1):
    T[i,ny] = 100.0*sin(pi*x[i]/w)

eps = 1e-8
itr = gauss_seidel(T,nx,ny,beta,eps)

# obtained solution in standard output
print "solution =",T
print "max iteration",itr
print "elapsed time =",time.time()-start

# plot setting
X, Y = meshgrid(x,y,indexing='ij')
contour(X,Y,T)
xlabel('X')
ylabel('Y')
cbar=colorbar()
cbar.set_label('T')
show()

