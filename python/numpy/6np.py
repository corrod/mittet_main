#coding utf-8

import numpy as np
from matplotlib import pyplot as plt

# broad cast

a = np.tile(np.arange(0, 40, 10), (3, 1)).T
print a
b = np.array([0, 1, 2])
print b
print a + b

d = np.arange(0, 40, 10)
print d.shape
print d
d = d[:, np.newaxis]
e = d[np.newaxis, :]
print d.shape
print d
print  b + d
