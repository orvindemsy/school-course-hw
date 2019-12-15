from numpy import *
from matplotlib.pyplot import *

# his program simpy draw the locus of amplification factor
# the deault conidtion when run locus will be inside the circle so it's stable
# if you change it to 1.0 locus will be outside unit cricle, so unstable, for FTCS method
# but if you change it to BTCS then it's stable, BTCS is always unstable


def FTCS(f, th, d):
    fm1 = f*exp(-1j*th)
    fp1 = f*exp(1j*th)
    # fp is plus fm is minus
    fnew = f + d*(fp1 - 2.0*f + fm1)
    return fnew


def BTCS(f, th, d):
    fm1 = f*exp(-1j*th)
    fp1 = f*exp(1j*th)
    fold = -d*fp1 + (1.0+2.0*d)*f - d*fm1
    fnew = (f/fold)*f
    return fnew


d = 0.5    # diffusion number
imax = 100
theta = linspace(0, 2*pi, imax)
x = zeros(imax, float)
y = zeros(imax, float)

xc = cos(theta)   # unit circle for reference
yc = sin(theta)

for i in range(0, imax):
    f = 1.0 + 0j    # arbitrary complex number input
    fnew = FTCS(f, theta[i], d)
    G = fnew/f
    x[i] = G.real
    y[i] = G.imag

# plot setting
xlabel('Re(G)')
ylabel('Im(G)')
plot(xc, yc, color='orange')
plot(x, y, marker='.', markersize=10, color='blue')
show()
