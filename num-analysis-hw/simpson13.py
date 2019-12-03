from numpy import *

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

a = 3.1
b = 3.9
imax = 9
x = zeros(imax, float)
f = zeros(imax, float)

for i in range(0, imax):
  x[i] = a + (b-a)/float(imax-1)*i
  f[i] = 1.0/x[i]
  print (i, ",", x[i], ",", f[i])

result = simpson13(x, f)
print("sum =", result)
err = abs(result-f[8])
print("error =", err)
