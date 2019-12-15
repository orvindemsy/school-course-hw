'''
Created by: Orvin Demsy
Student ID: B9TM9122

Exercise for 1st November lecture.
Please tackle the following projects when solving Laplace equation.
1. Make a plot for number of grids vs. elapsed time in a log-log scale by changing number of grids with three different
    methods: Gauss-Seidel, Gauss elimination, and conjugate gradient method.
2. Make a figure of contours by changing the right boundary value to x = 5 and the boundary condition at the right
    boundary to Neumann condition of dT/dx = 0 with one of the three methods.

Number 1
'''

import laplace_solver
from matplotlib import pyplot as plt

# Matrix containing time elapsed of each methods
time_GE = []
time_CG = []
time_GS = []

# Matrix containing number of grid
no_grid =[]

for i in (10, 20, 30):
    t1 = laplace_solver.laplace_GE(i, i)
    t2 = laplace_solver.laplace_CG(i, i)
    t3 = laplace_solver.laplace_GS(i, i)

    time_GE.append(t1)
    time_CG.append(t2)
    time_GS.append(t3)
    no_grid.append(i)

print(time_GE)
print(time_GS)
print(time_CG)
print(no_grid)
#
plt.figure()
plt.title('Log-Log Plot of Process Time of Each Method vs. Number of Grids')
plt.xlabel('Number of Grids')
plt.ylabel('Process Time of Each Method')
plt.loglog(no_grid, time_CG, basex=10, basey=10, marker='o', markersize='7')
plt.loglog(no_grid, time_GE, basex=10, basey=10, marker='x', markersize='7')
plt.loglog(no_grid, time_GS, basex=10, basey=10, marker='s', markersize='7')
plt.gca().legend(('Conjugate Gradient','Gauss Elimination','Gauss Seidel'))
plt.grid(True, which='both', ls=':')

plt.show()