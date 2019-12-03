from numpy import *
import matplotlib.pyplot as plt
from matplotlib import animation
# first mode 2 FTCS
# we can see the amplitude is growing, the FTCS is not unstable
#
# the amplitude is damped but solution stable, the orange one is the solution
# upwind is damped but better than lax
# but if we used upwind with higher mode
# set mode = 10, it will rapidly damped

# Lax wendroff scheme please make it
# describe in handout 11.37 page 666


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
plt.title('theta = %1.5f, c = %1.1f' % (theta, c))
plt.xlim(0, imax-1)
plt.xlabel('grid number')
plt.ylabel('f')
img = plt.plot(xr, fr, color='orange') \
    + plt.plot(x, f, marker='.', markersize=10, color='blue')
ims.append(img)

for itr in range(1, itrmax+1):
    f[0] = f[imax]      # periodic boundary condition
    f[imax+1] = f[1]    # periodic boundary condition
    f = FTCS(f, c, imax)
    fr = cos(theta*(xr-c*itr))   # exact solution
    img = plt.plot(xr, fr, color='orange')\
        + plt.plot(x, f, marker='.', markersize=10, color='blue')
    ims.append(img)

ani = animation.ArtistAnimation(fig, ims, interval=500)
plt.show()
