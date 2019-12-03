'''
Created by: Orvin Demsy
Student ID: B9TM9122

Exercise for November 15th lecture. Please make plots of von Neumann analysis of Lax-Wendroff scheme for convection
problem with convection numbers of 0.5, 1.0, and 1.5. In addition, make a plot to show time variation of numerical
solutions obtained by Lax-Wendroff scheme with arbitrary mode (wave number) and
convection number in the convection problem.

Task 2
'''


from numpy import *
import matplotlib.pyplot as plt
from matplotlib import animation

# FDE expression
def FTCS(f, c, imax):
    fnew = zeros(imax+2, float)
    for i in range(1, imax+1):
        fnew[i] = f[i] - c*0.5*(f[i+1] - f[i-1])
    return fnew


def lax(f, c, imax):
    fnew = zeros(imax+2, float)
    for i in range(1, imax+1):
        fnew[i] = 0.5*(f[i+1] + f[i-1]) - c*0.5*(f[i+1] - f[i-1])
    return fnew


def upwind(f, c, imax):
    fnew = zeros(imax+2, float)
    for i in range(1, imax+1):
        fnew[i] = f[i] - c*(f[i] - f[i-1])
    return fnew

def lax_wendroff(f, c, imax):
    fnew = zeros(imax + 2, float)
    for i in range(1, imax + 1):
        fnew[i] = f[i] - c * 0.5 * (f[i+1] - f[i-1]) + c**2 * 0.5 * (f[i+1] - 2*f[i] + f[i-1])
    return fnew

imax = 50
itrmax = 20

c = 0.5    # convection number
# mode number change to higher so larger number of wave
# signle mode mean same as previous, single frequency

m = 2      # mode number
theta = 2*pi*m/imax   # phase difference in single grid interval

x = zeros(imax+2, float)     # initial condition
f = zeros(imax+2, float)
for i in range(0, imax+2):
    x[i] = i-1
    f[i] = cos((i-1)*theta)

xr = linspace(0, imax, 1000)   # exact solution for reference
fr = cos(theta*xr)

fig = plt.figure()
ims = []
plt.title('Lax-wendroff Time Variation Plot \n theta = %1.5f, c = %1.1f' % (theta, c))
plt.xlim(0, imax-1)
plt.xlabel('grid number')
plt.ylabel('f')
img = plt.plot(xr, fr, color='orange') \
    + plt.plot(x, f, marker='.', markersize=10, color='blue')
ims.append(img)

for itr in range(1, itrmax+1):
    f[0] = f[imax]      # periodic boundary condition
    f[imax+1] = f[1]    # periodic boundary condition
    f = lax_wendroff(f, c, imax)
    fr = cos(theta*(xr-c*itr))   # exact solution
    img = plt.plot(xr, fr, color='orange')\
        + plt.plot(x, f, marker='.', markersize=10, color='blue')
    ims.append(img)

ani = animation.ArtistAnimation(fig, ims, interval=500)
plt.show()
