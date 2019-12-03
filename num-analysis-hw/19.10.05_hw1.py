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
    print("y({}) is {}".format(a, y(a)))
    print("y({}) is {}".format(b, y(b)))

    if y(a) > 0 and y(b) > 0 or y(a) < 0 and y(b) < 0:
        print("\nTry inputting different value\n")
    else:
        epsilon = float(input("Enter desired epsilon for accuracy: "))
        #Iterating function starts here
        iter = 1
        while (abs(a-b) > epsilon):
            c = (a + b)/2
            if y(a)*y(c) < 0:
                print("\n=========================================")
                print("This is iteration number %d" % iter)
                print("Value of y(a) or y(%.2f) is %.3f" % (a, y(a)))
                print("Value of y(b) or y(%.2f) is %.3f" % (b, y(b)))
                print("Value of y(c) or y(%.2f) is %.3f" % (c, y(c)))
                a = a
                b = c
                print("Our new a is %.3f" % a)
                print("Our new b is %.3f" % b)
                print("Current value of root/c is %.3f" % abs(c))
                print("=========================================")
                iter += 1
            elif y(c)*y(b) < 0:
                print("\n=========================================")
                print("This is iteration number %d" % iter)
                print("Value of y(a) or y(%.2f) is %.3f" % (a, y(a)))
                print("Value of y(b) or y(%.2f) is %.3f" % (b, y(b)))
                print("Value of y(c) or y(%.2f) is %.3f" % (c, y(c)))
                a = c
                b = b
                print("Our new a is %.3f" % a)
                print("Our new b is %.3f" % b)
                print("Current value of root/c is %.3f" % abs(c))
                print("=========================================")
                iter += 1
            else:
                print("Something's wrong check you initial a and b value")