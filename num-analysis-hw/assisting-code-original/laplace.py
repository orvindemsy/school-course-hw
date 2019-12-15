#!/usr/bin/env python3

# A sample code for numerically solving Laplace equation with 2nd-order centered-space approximation
# This code does not include convergence check in the Gauss-Seidel iteration, so you cannot obtain an accurate solution without any modification.
# For comparing with an analytic solution, you can add the exact values in the plot using mathematical function provided in the Python library.

from numpy import *
from matplotlib.pyplot import *

def gauss_seidel():
# definition of physical domain
    w = 10.0
    h = 15.0

# definition of computational domain
# i's are defined in [0,Nx], but i=0 and i=Nx are only used for boundary condition
    Nx = 10
# j's are defined in [0,Ny], but j=0 and j=Ny are only used for boundary condition
    Ny = 10

# maximum number of iteration (larger number is required for obtaining an accurate solution)
    itrmax = 10

# grid definition
    x = linspace(0, w, Nx+1)
    y = linspace(0, h, Ny+1)
    dx = w/Nx
    dy = h/Ny
    dx2i = 1.0/(dx*dx)
    dy2i = 1.0/(dy*dy)

# initialization of T
    T = zeros((Nx+1,Ny+1))

# given boundary condition for upper boundary (j = Ny)
    for i in range(0, Nx+1):
        T[i,Ny] = 100.0*sin(pi*x[i]/w)

# start Gauss-Seidel procedure
    for itr in range(0, itrmax):
        for j in range(1, Ny):
            for i in range(1, Nx):
                T[i,j]=((T[i-1,j]+T[i+1,j])*dx2i+(T[i,j-1]+T[i,j+1])*dy2i)/(2.0*dx2i+2.0*dy2i)

# obtained solution in standard output
    print(T)

# plot setting
    X, Y = meshgrid(x,y,indexing='ij')
    contour(X,Y,T)
    xlabel('X')
    ylabel('Y')
    cbar=colorbar()
    cbar.set_label('T')
    show()

if __name__ == '__main__':
    print('Solving Laplace equation with 2nd-order centered-space approximation...')
# Gauss-Seidel method
    gauss_seidel()
