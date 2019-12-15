'''
Created by: Orvin Demsy
Student ID: B9TM9122

Exercise for November 15th lecture. Please make plots of von Neumann analysis of Lax-Wendroff scheme for convection
problem with convection numbers of 0.5, 1.0, and 1.5. In addition, make a plot to show time variation of numerical
solutions obtained by Lax-Wendroff scheme with arbitrary mode (wave number) and
convection number in the convection problem.

Task 1
'''


from numpy import *
import matplotlib.pyplot as plt
from matplotlib import animation

def FTCS1(f, th, c):
   fm1 = f*exp(-1j*th)
   fp1 = f*exp( 1j*th)
   fnew = f - c*0.5*(fp1 - fm1)
   return fnew

def lax1(f, th, c):
   fm1 = f*exp(-1j*th)
   fp1 = f*exp( 1j*th)
   fnew = 0.5*(fp1 + fm1) - c*0.5*(fp1 - fm1)
   return fnew

def upwind1(f, th, c):
   fm1 = f*exp(-1j*th)
   fnew = f - c*(f - fm1)
   return fnew


def lax_wendroff1(f, th, c):
   fm1 = f * exp(-1j * th)
   fp1 = f * exp(1j * th)
   fnew = f - (c * 0.5) * (fp1 - fm1) + (c**2 * 0.5) * (fp1 - 2*f + fm1)
   return fnew

c = 0.5    # convection number
c = round(c, 2)
imax = 100
theta = linspace(0, 2*pi, imax)
x = zeros(imax, float)
y = zeros(imax, float)

xc = cos(theta) # unit circle for reference
yc = sin(theta)

for i in range(0, imax):
   f = 1.0 + 0j    # arbitrary complex number input
   fnew = lax_wendroff1(f, theta[i], c)
   G = fnew/f
   x[i] = G.real
   y[i] = G.imag

# plot setting
plt.figure('von Neumann Analysis of Lax-Wendroff')
plt.title('von Neumann Analysis of Lax-Wendroff \n convection number: %.2f' %c)
plt.xlabel('Re(G)')
plt.ylabel('Im(G)')
plt.plot(xc, yc, color='orange')
plt.plot(x, y, marker='.', markersize=10, color='blue')
plt.show()