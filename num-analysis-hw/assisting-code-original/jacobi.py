from numpy import *

def jacobi(a,b,eps,imax):
   n = b.size
   x = zeros(n, float)
   xnew = zeros(n, float)
   for it in range(0,imax):
      dxmax = 0.0
      for i in range(0,n):
         res = b[i]
         for j in range(0,n):
            res -= a[i,j]*x[j]
         dxmax = max(abs(res),dxmax)
         xnew[i] = x[i] + res/a[i,i]
      x = xnew
      print it, dxmax
      if dxmax < eps: return x

a = array([[ 4.0,-1.0, 0.0, 1.0, 0.0],
           [-1.0, 4.0,-1.0, 0.0, 1.0],
           [ 0.0,-1.0, 4.0,-1.0, 0.0],
           [ 1.0, 0.0,-1.0, 4.0,-1.0],
           [ 0.0, 1.0, 0.0,-1.0, 4.0]])
b = array([100.0,100.0,100.0,100.0,100.0])
imax = 100
eps = 1e-13

print a,",",b
sol = jacobi(a,b,eps,imax)
print sol

