from numpy import *

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

for i in range(0, imax):
  x[i] = a + (b-a)/float(imax-1)*i
  f[i] = 1.0/x[i]
  print(i, ",", x[i], ",", f[i])

result = trapezoid(x,f)
print("sum =", result)