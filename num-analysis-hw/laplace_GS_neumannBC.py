'''
Created by: Orvin Demsy
Student ID: B9TM9122

Exercise for 1st November lecture.
Please tackle the following projects when solving Laplace equation.
1. Make a plot for number of grids vs. elapsed time in a log-log scale by changing number of grids with three different
    methods: Gauss-Seidel, Gauss elimination, and conjugate gradient method.
2. Make a figure of contours by changing the right boundary value to x = 5 and the boundary condition at the right
    boundary to Neumann condition of dT/dx = 0 with one of the three methods.

Number 2
'''

from numpy import *
from matplotlib.pyplot import *
import time

# From eq 9.24
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

# Added time elapsed function
def time_elapsed(t0):
    t1 = time.time()
    return t1-t0

start = time.time()

# definition of computational domain
w = 5.0
h = 15.0
nx = 10 # grid for x direction
ny = 10 # grid for y direction
dx = w/nx
dy = h/ny
beta = dx/dy

# initialization
x = linspace(0, w, nx+1)
y = linspace(0, h, ny+1)
T = zeros((nx+1,ny+1), float)

# Dirichlet boundary condition for upper boundary (j = ny)
# Neumann boundary condition for rightside boundary (x = nx)
for i in range (0, ny+1):
    # We put the derivative of 100.0*sin(pi*x[i]/w) here
    T[nx, i] = (100.0*pi*cos((pi*x[i]/w)))/w

eps = 1e-8
itr = gauss_seidel(T,nx,ny,beta,eps)

# obtained solution in standard output
time_e = time_elapsed(start)
print ("solution =",T)
print ("max iteration",itr)
# print ("elapsed time =",time.time()-start)
print ("elapsed time =", time_e)

# plot setting
X, Y = meshgrid(x,y,indexing='ij')
contour(X,Y,T)
xlabel('X')
ylabel('Y')
cbar=colorbar()
cbar.set_label('T')
show()

