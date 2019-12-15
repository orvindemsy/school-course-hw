from numpy import *
from matplotlib.pyplot import *

# the default when run is outside the unit circle, so FTCS is unstable
# for convention 0.5
# FTCS is always unstable for convection
# for lax method is stable for 0.5
# fro 1.0 locus is coinside w unit circle, so still unstable
# but 1.5 is outside, so lax method is unstable
# for upwind 1.5 convec number, you can find smaller number inside circle, so smaller phase error, more diffusinion
# because ampitude is smaller compare to lax, upwind is more diffusing

# last he ran 1.0
# the locus coincide with unit circle
# with 0.5 is outside unit circle


def FTCS(f, th, c):
    fm1 = f*exp(-1j*th)
    fp1 = f*exp(1j*th)
    fnew = f - c*0.5*(fp1 - fm1)
    return fnew


def lax(f, th, c):
    fm1 = f*exp(-1j*th)
    fp1 = f*exp(1j*th)
    fnew = 0.5*(fp1 + fm1) - c*0.5*(fp1 - fm1)
    return fnew


def upwind(f, th, c):
    fm1 = f*exp(-1j*th)
    fnew = f - c*(f - fm1)
    return fnew


c = 0.5    # convection number
imax = 100
theta = linspace(0, 2*pi, imax)
x = zeros(imax, float)
y = zeros(imax, float)

xc = cos(theta)  # unit circle for reference
yc = sin(theta)

for i in range(0, imax):
    f = 1.0 + 0j    # arbitrary complex number input
    fnew = FTCS(f, theta[i], c)
    G = fnew/f
    x[i] = G.real
    y[i] = G.imag

# plot setting
xlabel('Re(G)')
ylabel('Im(G)')
plot(xc, yc, color='orange')
plot(x, y, marker='.', markersize=10, color='blue')
show()
