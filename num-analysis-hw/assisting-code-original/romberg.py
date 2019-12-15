from numpy import *

def trapezoid(x,f):
   n = x.size
   sum = f[0]+f[n-1]
   for i in range(1,n-1):
      sum += 2.0*f[i]
   sum = sum*(x[n-1]-x[0])/float(n-1)/2.0
   return sum

def romberg(n,x,f):
   imax = x.size
   s = zeros((n,n), float)
   for i in range(0,n):
      xt = zeros(2**i+1, float)
      ft = zeros(2**i+1, float)
      it = 0
      for j in range(0,imax,2**(n-1-i)):
         xt[it] = x[j]
         ft[it] = f[j]
         it += 1
      s[i,0] = trapezoid(xt,ft)
   for j in range(1,n):
      ex = float(2*j)
      den = 2.0**ex-1.0
      k = n-j
      for i in range(0,k):
         s[i,j] = s[i+1,j-1]+(s[i+1,j-1]-s[i,j-1])/den
   print s
   return s[0,n-1]

a = 3.1
b = 3.9
num = 4
imax = 2**(num-1)+1
x = zeros(imax, float)
f = zeros(imax, float)

for i in range(0,imax):
  x[i] = a + (b-a)/float(imax-1)*i
  f[i] = 1.0/x[i]
  print i,",",x[i],",",f[i]

result = romberg(num,x,f)
print "sum =",result

