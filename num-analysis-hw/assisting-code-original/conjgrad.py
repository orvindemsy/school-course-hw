from numpy import *

def conjgrad(a,b,eps,imax):
   n = b.size
   x = zeros(n, float)
   r0 = b - dot(a,x)
   p = r0
   for i in range(0, imax):
      alpha = dot(r0,r0)/dot(p,dot(a,p))
      x = x + alpha*p
      r1 = r0 - alpha*dot(a,p)
      res = linalg.norm(r1)
      beta = dot(r1,r1)/dot(r0,r0)
      p = r1 + beta*p
      r0 = r1
      print i, res
      if res < eps: return x

a = array([[ 4.0,-1.0, 0.0, 1.0, 0.0],
           [-1.0, 4.0,-1.0, 0.0, 1.0],
           [ 0.0,-1.0, 4.0,-1.0, 0.0],
           [ 1.0, 0.0,-1.0, 4.0,-1.0],
           [ 0.0, 1.0, 0.0,-1.0, 4.0]])
b = array([100.0,100.0,100.0,100.0,100.0])
imax = 100
eps = 1e-13

print a,",",b
sol = conjgrad(a,b,eps,imax)
print sol

