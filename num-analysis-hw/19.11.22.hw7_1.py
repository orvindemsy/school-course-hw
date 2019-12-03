'''
Created by: Orvin Demsy
Student ID: B9TM9122

Final assignment is a project for Laplace equation solver based
on finite element method (FEM). Please tackle at least one of the followings.
1. (Minimum) Discuss an order of accuracy of the FEM Laplace equation solver by making a log-log plot of numerical error
 as a function of element size. Definition of the numerical error and element size should be described in your document
 since there are some different definitions.
2. (Standard) Solve the same problem for Neumann boundary condition given in the second exercise of Nov. 1st and show
the nodal equation at the right boundary.
3. (Advanced) Upgrade the program code of FEM for advanced problems by employing nonuniform element size, triangular
element, higher-oder shape function, and so on. Show problem to be solved, result, and modification in the code.
'''

from numpy import *
from matplotlib import pyplot as plt


# The gauss_seidel function is obtain from the given laplace_FEM.py
def gauss_seidel(eps, imax, jmax):
   nmax = imax * jmax  # number of nodes including boundary
   emax = (imax - 1) * (jmax - 1)  # number of elements
   coef = zeros((nmax, nmax), float)  # coefficinets of nodal equations
   rhs = zeros(nmax, float)  # RHS of nodal equations
   ele = zeros((emax, 4), integer)  # set of node indices of element

   # problem setting
   xmax = 10.0
   xmin = 0.0
   ymax = 15.0
   ymin = 0.0
   dx = (xmax - xmin) / (imax - 1)
   dy = (ymax - ymin) / (jmax - 1)
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
         x[na + nb] = xmin + i * dx
         y[na + nb] = ymin + j * dy
         if (i == 0 or j == 0 or i == imax - 1 or j == jmax - 1):  # Dirichlet boundary
            if (j == jmax - 1):
               T[na + nb] = 100.0 * sin(pi * x[na + nb] / (xmax - xmin))
            else:
               T[na + nb] = 0.0
            ib.append(na + nb)  # indices of boundary nodes
            nb += 1  # number of boundary nodes
         else:
            ia.append(na + nb)  # indices of active nodes
            na += 1  # number of active nodes

   # definition of elements as a set of nodes
   ne = 0
   for j in range(0, jmax - 1):
      for i in range(0, imax - 1):
         ele[ne, 0] = i + j * imax  # anti-clockwise order
         ele[ne, 1] = i + 1 + j * imax
         ele[ne, 2] = i + 1 + (j + 1) * imax
         ele[ne, 3] = i + (j + 1) * imax
         ne += 1

   # creating matrix coefficient for nodal equation from element equations
   for ne in range(0, emax):
      n0 = ele[ne, 0]
      n1 = ele[ne, 1]
      n2 = ele[ne, 2]
      n3 = ele[ne, 3]
      # W1=N1
      c0 = -1.0 / 3.0 * (dy / dx + dx / dy)
      c1 = 1.0 / 6.0 * (2.0 * dy / dx - dx / dy)
      c2 = 1.0 / 6.0 * (dy / dx + dx / dy)
      c3 = 1.0 / 6.0 * (-dy / dx + 2.0 * dx / dy)
      coef[n0, n0] += c0
      coef[n0, n1] += c1
      coef[n0, n2] += c2
      coef[n0, n3] += c3
      if (n1 in ib): rhs[n0] += -c1 * T[n1]  # known boundary values to be in RHS
      if (n2 in ib): rhs[n0] += -c2 * T[n2]
      if (n3 in ib): rhs[n0] += -c3 * T[n3]
      # W2=N2
      c0 = 1.0 / 6.0 * (2.0 * dy / dx - dx / dy)
      c1 = -1.0 / 3.0 * (dy / dx + dx / dy)
      c2 = 1.0 / 6.0 * (-dy / dx + 2.0 * dx / dy)
      c3 = 1.0 / 6.0 * (dy / dx + dx / dy)
      coef[n1, n0] += c0
      coef[n1, n1] += c1
      coef[n1, n2] += c2
      coef[n1, n3] += c3
      if (n0 in ib): rhs[n1] += -c0 * T[n0]
      if (n2 in ib): rhs[n1] += -c2 * T[n2]
      if (n3 in ib): rhs[n1] += -c3 * T[n3]
      # W3=N3
      c0 = 1.0 / 6.0 * (dy / dx + dx / dy)
      c1 = 1.0 / 6.0 * (-dy / dx + 2.0 * dx / dy)
      c2 = -1.0 / 3.0 * (dy / dx + dx / dy)
      c3 = 1.0 / 6.0 * (2.0 * dy / dx - dx / dy)
      coef[n2, n0] += c0
      coef[n2, n1] += c1
      coef[n2, n2] += c2
      coef[n2, n3] += c3
      if (n0 in ib): rhs[n2] += -c0 * T[n0]
      if (n1 in ib): rhs[n2] += -c1 * T[n1]
      if (n3 in ib): rhs[n2] += -c3 * T[n3]
      # W4=N4
      c0 = 1.0 / 6.0 * (-dy / dx + 2.0 * dx / dy)
      c1 = 1.0 / 6.0 * (dy / dx + dx / dy)
      c2 = 1.0 / 6.0 * (2.0 * dy / dx - dx / dy)
      c3 = -1.0 / 3.0 * (dy / dx + dx / dy)
      coef[n3, n0] += c0
      coef[n3, n1] += c1
      coef[n3, n2] += c2
      coef[n3, n3] += c3
      if (n0 in ib): rhs[n3] += -c0 * T[n0]
      if (n1 in ib): rhs[n3] += -c1 * T[n1]
      if (n2 in ib): rhs[n3] += -c2 * T[n2]

   # deleting boundary nodes from matrix elements and RHS vector
   A = delete(coef, ib, 0)
   A = delete(A, ib, 1)
   b = delete(rhs, ib, 0)

   itrmax = 1000
   f = zeros(b.size, float)   # initial guess
   for itr in range(0, itrmax):
      err = 0.0
      for n in range(0, b.size):
         new = (-dot(A[n,:],f)+A[n,n]*f[n]+b[n])/A[n,n]
         err = max(abs(f[n]-new),err)
         f[n] = new
      if (err < eps): return f, itr, A, b

   return f, itrmax, A, b

# Generate approximate value with grid 5 x 7,
imax1 = 5 # including boundary
jmax1 = 7 # including boundary
eps = 1e-8
f1, itr1, A1, b1 = gauss_seidel(eps, imax1, jmax1)

# Print the approximate result
# print("\nThe approximation solution of Laplace Equation by FEM with grid 5 x 7: \n", f1)

# Matching the values of f with table 12.7 in page 751,
# it is inferred that at x = 2.5 and y = 12.5 the approximate value is 30.8511

# The error at (2.5, 12.5) and (5.0, 12.5) in regard to grid will be observed.
exact1 = 32.229780
exact2 = 45.579791
print("The exact value at (2.5, 12.5) is :", exact1)
print("The exact value at (5.0, 12.5) is :", exact1)

# The error at (2.5, 12.5) with grid 5 x 7 is
print("\nThe approx value at (2.5, 12.5) with grid 5 x 7 is: ", f1[12])
print("The error at (2.5, 12.5) is", abs(exact1-f1[12]))

# The error at (5.0, 12.5) with grid 5 x 7 is
print("\nThe approx value at (5.0, 12.5) with grid 5 x 7 is: ", f1[13])
print("The error at (5.0, 12.5) is", abs(exact2-f1[13]))

# ======================================================================
# Generate approximate value with grid 9 x 13
imax2 = 9 # including boundary
jmax2 = 13 # including boundary
f2, itr2, A2, b2 = gauss_seidel(eps, imax2, jmax2)

# Print the approximate result
# print("\nThe approximation solution of Laplace Equation by FEM with grid 9 x 13: \n", f2)

# The error at (2.5, 12.5) with grid 9 x 13 is
print("\nThe approx value at (2.5, 12.5) with grid 9 x 13 is: ", f2[64])
print("The error at (2.5, 12.5) is", abs(exact1-f2[64]))

# The error at (5.0, 12.5) with grid 9 x 13 is
print("\nThe approx value at (5.0, 12.5) with grid 9 x 13 is: ", f2[66])
print("The error at (5.0, 12.5) is", abs(exact2-f2[66]))
# ======================================================================
# Generate approximate value with grid 17 x 25
imax3 = 17 # including boundary
jmax3 = 25 # including boundary
f3, itr3, A3, b3 = gauss_seidel(eps, imax3, jmax3)

# Print the approximate result
# print("\nThe approximation solution of Laplace Equation by FEM with grid 17 x 25: \n", f3)

# The error at (2.5, 12.5) with grid 9 x 13 is
print("\nThe approx value at (2.5, 12.5) with grid 17 x 25 is: ", f3[288])
print("The error at (2.5, 12.5) is", abs(exact1-f3[288]))

# The error at (5.0, 12.5) with grid 9 x 13 is
print("\nThe approx value at (5.0, 12.5) with grid 17 x 25 is: ", f3[292])
print("The error at (5.0, 12.5) is", abs(exact2-f3[292]))
# ======================================================================
# Creating plot
# For (2.5, 12.5)
num_error1 = [abs(exact1-f1[12]), abs(exact1-f2[64]), abs(exact1-f3[288])]
elem_size1 = [(imax1*jmax1), (imax2*jmax2), (imax3*jmax3)]

# For (5.0, 12.5)
num_error2 = [abs(exact2-f1[13]), abs(exact2-f2[66]), abs(exact2-f3[292])]
elem_size2 = [(imax1*jmax1), (imax2*jmax2), (imax3*jmax3)]
#
plt.figure()
plt.title('Log-Log Plot of \n Number of Grids (Element Size) vs Error of FEM Laplace Solution')
plt.xlabel('Number of Grids (Element Size) ')
plt.ylabel('Error')
plt.loglog(elem_size1, num_error1, basex=10, basey=10, marker='o', markersize='7')
plt.loglog(elem_size2, num_error2, basex=10, basey=10, marker='o', markersize='7')
# plt.loglog(no_grid, time_GE, basex=10, basey=10, marker='x', markersize='7')
plt.gca().legend(('at x = 2.5 y = 12.5','at x = 5.0 y = 12.5'))
plt.grid(True, which='both', ls=':')

plt.show()
