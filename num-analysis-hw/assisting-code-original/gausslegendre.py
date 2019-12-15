from numpy import *

def f(x,n):
   p = zeros(n, float)
   p1 = 1.0
   p2 = 0.0
   for i in range(0,n):
      p3 = p2
      p2 = p1
      p1 = ((2.0*float(i+1)-1.0)*x*p2-float(i)*p3)/float(i+1)
      p[i] = p1
   return p

def gausslegendre(n):
   z = zeros(n, float)
   w = zeros(n, float)
   m = int((n+1)/2)
   eps = 1e-8
   for i in range(0,m):
      x = cos(pi*(float(i)+0.75)/(float(n)+0.5))
      for it in range(0,100):
         p = f(x,n)
         pp = float(n)*(x*p[n-1]-p[n-2])/(x*x-1.0)
         x1 = x
         x = x1-p[n-1]/pp
         if abs(x-x1) < eps: break
      z[i] = -x
      z[n-1-i] = x
      w[i] = 2.0/((1.0-z[i]*z[i])*pp*pp)
      w[n-1-i]=w[i]
   return z,w

n = 2
z,w = gausslegendre(n)
print "Abscissas =",z
print "Weights =",w

a = 3.1
b = 3.9
c = 0.5*(b+a)
m = 0.5*(b-a)
z = c - m*z
w = w*m
sum = 0.0
for i in range(0,n):
   sum += w[i]*(1.0/z[i])
print "Integration of (1/x) =",sum

