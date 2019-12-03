from numpy import *
from matplotlib.pyplot import *

def gauss_seidel(A, b, eps):
   itrmax = 1000
   f = zeros(b.size, float)   # initial guess
   for itr in range(0, itrmax):
      err = 0.0
      for n in range(0, b.size):
         new = (-dot(A[n,:],f)+A[n,n]*f[n]+b[n])/A[n,n]
         err = max(abs(f[n]-new),err)
         f[n] = new
      if (err < eps): return f, itr
   return f, itrmax

imax = 5    # including boundary
jmax = 7    # including boundary
nmax = imax*jmax   # number of nodes including boundary
emax = (imax-1)*(jmax-1)   # number of elements
coef = zeros((nmax,nmax), float)   # coefficinets of nodal equations
rhs = zeros(nmax, float)           # RHS of nodal equations
ele = zeros((emax,4), integer)     # set of node indices of element

# problem setting
xmax = 10.0
xmin = 0.0
ymax = 15.0
ymin = 0.0
dx = (xmax-xmin)/(imax-1)
dy = (ymax-ymin)/(jmax-1)
x = zeros(nmax, float)
y = zeros(nmax, float)
T = zeros(nmax, float)

# definition of nodes
na = 0
nb = 0
ia = []
ib = []
for j in range(0, jmax):
   for i in range(0, imax):
      x[na+nb] = xmin + i*dx
      y[na+nb] = ymin + j*dy
      if ( i==0 or j==0 or i==imax-1 or j==jmax-1 ): # Dirichlet boundary
         if (j==jmax-1):
            T[na+nb] = 100.0*sin(pi*x[na+nb]/(xmax-xmin))
         else:
            T[na+nb] = 0.0
         ib.append(na+nb)     # indices of boundary nodes
         nb += 1              # number of boundary nodes
      else:
         ia.append(na+nb)     # indices of active nodes
         na += 1              # number of active nodes

# definition of elements as a set of nodes
ne = 0
for j in range(0, jmax-1):
   for i in range(0, imax-1):
      ele[ne,0] = i+j*imax         # anti-clockwise order
      ele[ne,1] = i+1+j*imax
      ele[ne,2] = i+1+(j+1)*imax
      ele[ne,3] = i+(j+1)*imax
      ne += 1

# creating matrix coefficient for nodal equation from element equations
for ne in range(0, emax):
   n0 = ele[ne,0]
   n1 = ele[ne,1]
   n2 = ele[ne,2]
   n3 = ele[ne,3]
# W1=N1
   c0 =-1.0/3.0*(dy/dx+dx/dy)
   c1 = 1.0/6.0*(2.0*dy/dx-dx/dy)
   c2 = 1.0/6.0*(dy/dx+dx/dy)
   c3 = 1.0/6.0*(-dy/dx+2.0*dx/dy)
   coef[n0,n0] += c0
   coef[n0,n1] += c1
   coef[n0,n2] += c2
   coef[n0,n3] += c3
   if (n1 in ib): rhs[n0] +=-c1*T[n1]   # known boundary values to be in RHS
   if (n2 in ib): rhs[n0] +=-c2*T[n2]
   if (n3 in ib): rhs[n0] +=-c3*T[n3]
# W2=N2
   c0 = 1.0/6.0*(2.0*dy/dx-dx/dy)
   c1 =-1.0/3.0*(dy/dx+dx/dy)
   c2 = 1.0/6.0*(-dy/dx+2.0*dx/dy)
   c3 = 1.0/6.0*(dy/dx+dx/dy)
   coef[n1,n0] += c0
   coef[n1,n1] += c1
   coef[n1,n2] += c2
   coef[n1,n3] += c3
   if (n0 in ib): rhs[n1] +=-c0*T[n0]
   if (n2 in ib): rhs[n1] +=-c2*T[n2]
   if (n3 in ib): rhs[n1] +=-c3*T[n3]
# W3=N3
   c0 = 1.0/6.0*(dy/dx+dx/dy)
   c1 = 1.0/6.0*(-dy/dx+2.0*dx/dy)
   c2 =-1.0/3.0*(dy/dx+dx/dy)
   c3 = 1.0/6.0*(2.0*dy/dx-dx/dy)
   coef[n2,n0] += c0
   coef[n2,n1] += c1
   coef[n2,n2] += c2
   coef[n2,n3] += c3
   if (n0 in ib): rhs[n2] +=-c0*T[n0]
   if (n1 in ib): rhs[n2] +=-c1*T[n1]
   if (n3 in ib): rhs[n2] +=-c3*T[n3]
# W4=N4
   c0 = 1.0/6.0*(-dy/dx+2.0*dx/dy)
   c1 = 1.0/6.0*(dy/dx+dx/dy)
   c2 = 1.0/6.0*(2.0*dy/dx-dx/dy)
   c3 =-1.0/3.0*(dy/dx+dx/dy)
   coef[n3,n0] += c0
   coef[n3,n1] += c1
   coef[n3,n2] += c2
   coef[n3,n3] += c3
   if (n0 in ib): rhs[n3] +=-c0*T[n0]
   if (n1 in ib): rhs[n3] +=-c1*T[n1]
   if (n2 in ib): rhs[n3] +=-c2*T[n2]

# deleting boundary nodes from matrix elements and RHS vector
A = delete(coef,ib,0)
A = delete(A,ib,1)
b = delete(rhs,ib,0)

eps = 1e-8
f, itr = gauss_seidel(A, b, eps)
print (A,b,itr)
print(f)
for n in range(0, na):   # restoring to active nodes
   T[ia[n]] = f[n]

# plot setting
contour(x.reshape(jmax,imax),y.reshape(jmax,imax),T.reshape(jmax,imax))
xlabel('X')
ylabel('Y')
cbar=colorbar()
cbar.set_label('T')
show()

