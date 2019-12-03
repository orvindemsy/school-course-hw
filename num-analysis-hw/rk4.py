""" From "COMPUTATIONAL PHYSICS" & "COMPUTER PROBLEMS in PHYSICS"
    by RH Landau, MJ Paez, and CC Bordeianu (deceased)
    Copyright R Landau, Oregon State Unv, MJ Paez, Univ Antioquia, 
    C Bordeianu, Univ Bucharest, 2017. 
    Please respect copyright & acknowledge our work."""

# rk4.py 4th order Runge Kutta application wi built in rk4
		 
from numpy import *
import matplotlib.pylab as plt

#   Initialization
a = 0.
b = 10.
n = 1000
ydumb = zeros((2), float);    y = zeros((2), float)
fReturn = zeros((2), float);  k1 = zeros((2), float)
k2 = zeros((2), float);       k3 = zeros((2), float) 
k4 = zeros((2), float)
y[0] = 3.;   y[1] = -5.
t = a;       h = (b-a)/n;
tt = zeros(n+1, float); yy = zeros(n+1, float); vy = zeros(n+1, float)

def f(t,x):                       # Force function 
    fReturn[0] = x[1]
    fReturn[1] = -100.*x[0] -2.*x[1] + 10.*sin(3.*t)
    return fReturn

def rk4(t,h,N,y):
    k1 = [0]*(N)
    k2 = [0]*(N)
    k3 = [0]*(N)
    k4 = [0]*(N)
    fR = [0]*(N)
    ydumb = [0]*(N)
    fR = f(t, y)                     # Returns RHS's
    for i in range(0, N):
       k1[i] = h*fR[i]
    for i in range(0, N):
        ydumb[i] = y[i] + k1[i]/2.
    k2 = h*f(t+h/2., ydumb)
    for i in range(0, N):
        ydumb[i] = y[i] + k2[i]/2.
    k3 = h*f(t+h/2., ydumb)
    for i in range(0, N):
        ydumb[i] = y[i] + k3[i]
    k4 = h*f(t+h, ydumb)
    for i in range(0, N):
        y[i] = y[i] + (k1[i] + 2.*(k2[i] + k3[i]) + k4[i])/6.
    return y

ye = zeros(2, float)
ye[0] = 3.;   ye[1] = -5.
yye = zeros(n+1, float)
vye = zeros(n+1, float)

def euler(t,h,N,ye):
    fR = f(t,ye)
    for i in range(0, N):
        ye[i] = ye[i] + h*fR[i]
    return ye

i = 0
while (t < b):                         # Time loop
    if ((t + h) > b):
        h = b - t                      # Last step
    y = rk4(t,h,2,y)
    tt[i] = t
    yy[i] = y[0]
    vy[i] = y[1]
    ye = euler(t,h,2,ye)
    yye[i] = ye[0]
    vye[i] = ye[1]
    t = t + h
    i = i + 1

fig, axes = plt.subplots(nrows=1,ncols=2,figsize=(12,5))
axes[0].plot(tt,yy,label='rk4')
axes[0].plot(tt,yye,label='Euler')
axes[0].grid()
axes[0].set_xlabel('t')
axes[0].set_ylabel('x(t)')
axes[1].plot(yy,vy,label='rk4')
axes[1].plot(yye,vye,label='Euler')
axes[1].grid()
axes[1].set_xlabel('x(t)')
axes[1].set_ylabel('v(t)')
plt.legend()
plt.show()

