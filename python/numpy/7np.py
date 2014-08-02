#coding utf-8

import numpy as np
from matplotlib import pyplot as plt

a = np.array([[1, 2, 3], [4, 5, 6]])
print a
print a.ravel()
print a.T
print a.T.ravel()
print a.shape
b = a.ravel()
print b
print b.reshape((2, 3))

k = np.arange(100)

print k.reshape((20, 5))

print np.arange(15).reshape(3, 5).T
