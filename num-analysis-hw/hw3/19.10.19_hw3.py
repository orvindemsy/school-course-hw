'''
Created by: Orvin Demsy
Student ID: B9TM9122

Exercise for October 19th lecture. Please submit a figure plotting errors in numerical integrations of the trapezoid rule and
Simpson's 1/3 rule with different increments (different number of discrete points) and confirm the slope in log-log plots,
which corresponds to the order of accuracy determined by the trancation error.
'''

from numpy import *
import sympy as sym
from matplotlib import pyplot as plt


def simpson13(imax):
    # Define the list
    x = zeros(imax, float)
    f = zeros(imax, float)

    for i in range(0, imax):
        x[i] = a + (b - a) / float(imax - 1) * i
        f[i] = 1.0 / x[i]

    # Calculating the approximation
    n = x.size
    sum = f[0]+f[n-1]
    sum2 = 0.0
    sum4 = 0.0
    for i in range(2, n-1,2):
      sum2 += f[i]
    for i in range(1, n-1,2):
      sum4 += f[i]
    sum = (sum+2.0*sum2+4.0*sum4)*(x[n-1]-x[0])/float(n-1)/3.0
    return sum


def trapezoid(imax):
    # Define the list
    x = zeros(imax, float)
    f = zeros(imax, float)

    for i in range(0, imax):
        x[i] = a + (b - a) / float(imax - 1) * i
        f[i] = 1.0 / x[i]

    # Calculating the approximation
    n = x.size
    sum = f[0]+f[n-1]
    for i in range(1, n-1):
        sum += 2.0*f[i]
    sum = sum*(x[n-1]-x[0])/float(n-1)/2.0
    return sum


def actual_value(a, b):
    # Compute the exact integration value
    x = sym.symbols('x')
    res = sym.integrate(1/x, (x, a, b))
    return res


# Boundary of integration
a = 3.1
b = 3.9

# Creating list of trapezoid result with increment of 4, 2, and 1
list_err_trap = []
list_err_sim = []
list_res_trap = []
list_res_sim = []
incre_trap = []
incre_sim = []

# List all necessary value for trapezoid
for i in (5, 3, 2):
    trap_res = trapezoid(i)
    ex_res = actual_value(a, b)
    inc = abs((b-a))/(i-1)
    inc = round(inc, 2)

    list_res_trap.append(trap_res)
    list_err_trap.append(trap_res - ex_res)
    incre_trap.append(inc)

# List all necessary value for simpson
# Simpson 1/3 should have even n, therefore n = 1 can't be included
for i in (5, 3):
    sim_res = simpson13(i)
    ex_res = actual_value(a, b)
    inc = abs((b - a)) / (i - 1)
    inc = round(inc, 2)

    list_res_sim.append(sim_res)
    list_err_sim.append(sim_res - ex_res)
    incre_sim.append(inc)

# Drawing plot
plt.figure(1)
plt.xlabel('Increment')
plt.ylabel('Error')
plt.title('Figures of Trapezoid and Simpson''s 1/3 Error Result \n in Regard to Increment Used ')

plt.loglog(incre_trap, list_err_trap, basex=10, basey=10, color='r', marker='o', markersize=7)
plt.loglog(incre_sim, list_err_sim, basex=10, basey=10, color='b', marker="^", markersize=7)
plt.gca().legend(('Trapezoid','Simpson 1/3'))
plt.grid(True, which="both", ls='--')
plt.show()