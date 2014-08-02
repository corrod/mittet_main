#coding utf-8

import numpy as np
from matplotlib import pyplot as plt

# a = np.array([1, 2, 3]) + 1.5
# print a

# a = np.array([1., 2, 3])
# print a.dtype
# a[0] = 2.1
# print a

# p = np.poly1d([3, 2, -1])
# print p
x = np.linspace(0, 1, 20)
y = np.cos(x) + 0.3*np.random.rand(20)
p = np.poly1d(np.polyfit(x, y, 3))
t = np.linspace(0, 1, 200)
plt.plot(x, y, 'o', t, p(t), '-')
plt.show()
