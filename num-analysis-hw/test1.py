from matplotlib import pyplot as plt
import numpy as np
import random

w = 4
h = 3
d = 70
plt.figure(figsize=(w, h), dpi = d)

list1 = np.random.randint(10, size=10)
list2 = 0.2 * np.random.randint(8, size=10)
list3 = []
n = []
act_val = 0.8
list4 = []

# def insert(x, index, item):
#     x[:]

def push(x, y):
    x = x + [y]
    return x

# push(list4, 8)
# print(list4[:])


def loh(i):
    res = 4**2/(i+1)
    return res

for i in range(10):
    list3.append(loh(i))
    list4.append(i)

print(list3)
print(list4)

# print(list1)
# print(list2)

# x = [20, 200, 2000, 20000, 200000]
# y = [30, 300, 3000, 30000, 300000]

# plt.loglog(list3, list4, basex=10, basey=10)
# plt.show()
