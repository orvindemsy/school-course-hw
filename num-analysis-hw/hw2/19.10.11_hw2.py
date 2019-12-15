'''
Homework for Numerical Analysis Class
Submit by: 19-Oct-2019
Created by: Orvin Demsy

Create algorithm for SOR (Succesive Over Relaxation) method
This SOR is based on the Jacobi Iteration Method

'''

from numpy import *

#Through trial and error, it is found that giving weight of 1.1 sppeds up the process
w = 1.1

def jacobi(a,b,eps,imax, w):
   n = b.size
   x = zeros(n, float)
   xnew = zeros(n, float)
   f = open("answer_hw2.txt", 'w')
   f.write("The iteration\n")
   for it in range(0,imax):
       dxmax = 0.0
       for i in range(0,n):
           res = b[i]
           for j in range(0,n):
               res -= a[i,j]*x[j]
           dxmax = max(abs(res),dxmax)
           xnew[i] = x[i] + w * res/a[i,i]
       x = xnew
       f.write("number of iter %d %.15f\n" %(it, dxmax))
      # print (it, dxmax)
       if dxmax < eps:
           return x



a = array([[ 4.0,-1.0, 0.0, 1.0, 0.0],
           [-1.0, 4.0,-1.0, 0.0, 1.0],
           [ 0.0,-1.0, 4.0,-1.0, 0.0],
           [ 1.0, 0.0,-1.0, 4.0,-1.0],
           [ 0.0, 1.0, 0.0,-1.0, 4.0]])
b = array([100.0,100.0,100.0,100.0,100.0])
imax = 100
eps = 1e-13



sol = jacobi(a, b, eps, imax, w)
sol.tolist()
# print(sol)

f = open("answer_hw2.txt", 'a')
f.write("\nThe matrix to be iterated\n")
f.write('a: ' + str(a) + '\n\n' + 'b: ' + str(b) + '\n\n')
f.write("The answer of x matrix:\n[")

for t in sol:
    f.write("{:.3f}  ".format(t))
f.write("]")
f.close()