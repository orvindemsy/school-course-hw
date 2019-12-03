""" From "COMPUTATIONAL PHYSICS" & "COMPUTER PROBLEMS in PHYSICS"
    by RH Landau, MJ Paez, and CC Bordeianu (deceased)
    Copyright R Landau, Oregon State Unv, MJ Paez, Univ Antioquia, 
    C Bordeianu, Univ Bucharest, 2017. 
    Please respect copyright & acknowledge our work."""
    
# rk4Duffing.py solve ODE for Duffing Osc via rk4 & Matplotlib

import numpy as np, matplotlib.pylab as plt
from math import *

tt = np.zeros((2000),float);  yy = np.zeros((2000),float);
vy = np.zeros((2000),float);  y = [0]*(2); rhs = [0]*(2)
a = 0.5;  b = -0.5;  g = 0.02;  F = 0.0008;   w = 1.  # Initial values
# a = 0.0;  b = 1.0;  g = 0.04;  F = 0.2;   w = 1.      # SOLUTION for A, three solution for a forced Duffing oscillator
# a = 0.0;  b = 5.3;  g = 0.009;  F = 0.4;   w = 1.       # SOLUTION for B, Ueda oscillator #REFENRECE FROM LI
h = 0.1                                     # Time step
i = 0; y[0] = 0.09;  y[1] =  0         # Initial x, velocity

def f(t,y):                            # RHS function
    rhs[0] = y[1]
    rhs[1] = -2*g*y[1] - a*y[0] - b*y[0]**3 + F*cos(w*t)
    return rhs

def rk4(t, h, N, y):
    k1 = [0]*(N); k2 = [0]*(N);  k3 = [0]*(N);  k4 = [0]*(N)
    fvector = [0]*(N)
    ydumb = [0]*(N)
    fvector = f(t, y)                     # Returns RHS's
    for i in range(0, N):  k1[i] = h*fvector[i]
    for i in range(0, N):  ydumb[i] = y[i] + k1[i]/2.
    fvector = f(t+h/2., ydumb)
    for i in range(0, N):  k2[i] = h*fvector[i]
    for i in range(0, N):  ydumb[i] = y[i] + k2[i]/2.
    fvector = f(t+h/2., ydumb)
    for i in range(0, N):  k3[i] = h*fvector[i]
    for i in range(0, N):  ydumb[i] = y[i] + k3[i] 
    fvector = f(t+h, ydumb)
    for i in range(0, N):  k4[i] = h*fvector[i]
    for i in range(0, N):  y[i] = y[i] + (k1[i]+2.*(k2[i]+k3[i])+k4[i])/6.
    return y

for t in np.arange(0,200,h):                      # Time Loop
    tt[i] = t
    y = rk4(t,h,2,y)                         # Call rk4
    yy[i] = y[0]                                         # x(t)
    vy[i] = y[1]                                         # v(t)
    i = i+1

fig, axes = plt.subplots(nrows=1,ncols=2,figsize=(12,5))
axes[0].plot(tt[1000:],yy[1000:])    # 1000 to avoid transients
axes[0].grid()                                            # x(t)
axes[0].set_title('Duffing Oscillator x(t)')
axes[0].set_xlabel('t')
axes[0].set_ylabel('x(t)')
axes[1].plot(yy[1000:],vy[1000:])
axes[1].grid()
axes[1].set_title('Phase Space Orbits for Duffing Oscillator')
axes[1].set_xlabel('x(t)')
axes[1].set_ylabel('v(t)')
plt.show()

