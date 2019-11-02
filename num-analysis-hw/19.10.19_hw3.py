'''
Exercise for October 19th lecture. Please submit a figure plotting errors in numerical integrations of the trapezoid rule and
Simpson's 1/3 rule with different increments (different number of discrete points) and confirm the slope in log-log plots,
which corresponds to the order of accuracy determined by the trancation error.
'''

from numpy import *
import sympy as sym
from matplotlib import pyplot as plt

def simpson13(x,f):
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

def trapezoid(x,f):
   n = x.size
   sum = f[0]+f[n-1]
   for i in range(1, n-1):
      sum += 2.0*f[i]
   sum = sum*(x[n-1]-x[0])/float(n-1)/2.0
   return sum

a = 3.1
b = 3.9
imax = 9
x = zeros(imax, float)
f = zeros(imax, float)

print("Actual Value")
for i in range(0, imax):
  x[i] = a + (b-a)/float(imax-1)*i
  f[i] = 1.0/x[i]
  print(i, ",", x[i], ",", f[i])

def actual_value(a, b):
    x = sym.symbols('x')
    res = sym.integrate(1/x, (x, 3.1, 3.9))
    return res

result_trap = trapezoid(x, f)
print("sum =", result_trap)

result_simp = simpson13(x, f)
print("sum =", result_simp)

result_exact = actual_value(a, b)
print(actual_value(3.1, 3.9))

err = abs(result_trap-f[8])
print("error =", err)

# for i in range (0, imax):