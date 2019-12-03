'''
Homework for Numerical Analysis Class
Submit by: 12-Oct-2019
Created by: Orvin Demsy

Given an equation find the root using one of four following approximation approaches:
1. Interval halving (bisection)
2. False position (regula falsi)
3. Secant method
4. Newton's method
'''
#Interval halving (bisection) is chosen to solve this equation
import math

#Equation
def y(x):
    res = 5/3 * math.cos((40/180 * math.pi)) - 5/2 * math.cos((x/180 * math.pi)) + 11/6 - math.cos((40/180 * math.pi) - (x/180 * math.pi))
    return res

#Checking the value of boundary which are y(a) and y(b)
#They can't be both positive or negative
a = 1
b = 1
while (y(a) > 0 and y(b) > 0 or y(a) < 0 and y(b) < 0 or a == b):
    print("Please input the value of a and b with a < b")
    print("otherwise a result will be wrong\n")
    a = int(input("Input a (1st boundary): "))
    b = int(input("Input b (2nd boundary): "))
    print("y({}) is {:.6f}".format(a, y(a)))
    print("y({}) is {:.6f}".format(b, y(b)))

    if y(a) > 0 and y(b) > 0 or y(a) < 0 and y(b) < 0:
        print("\nTry inputting different value\n")
    else:
        epsilon = float(input("Enter desired epsilon for accuracy: "))
        #Iterating function starts here
        iter = 1
        print("Iter num\t\tValue of a\t\tValue of y(a)\t\tValue of b\t\tValue of y(b)\t\tValue of c\t\tValue of y(c)")
        while (abs(a-b) > epsilon):

            c = (a + b)/2
            if y(a)*y(c) < 0:
                print("%d\t\t\t\t%.6f\t\t%.6f\t\t\t%.6f\t\t%.6f\t\t\t%.6f\t\t%.6f" %(iter, a, y(a), b, y(b), c, y(c)))
                a = a
                b = c 
                iter += 1
            elif y(c)*y(b) < 0:
                print("%d\t\t\t\t%.6f\t\t%.6f\t\t\t%.6f\t\t%.6f\t\t\t%.6f\t\t%.6f" %(iter, a, y(a), b, y(b), c, y(c)))
                a = c
                b = b
                iter += 1
            else:
                print("Something's wrong check you initial a and b value")