#coding utf-8

import numpy as np
from matplotlib import pyplot as plt

from scipy import linalg
arr = np.array([[1., 2],
                [3, 4]])
iarr = linalg.inv(arr)
print arr
print iarr
print linalg.det(arr)
print np.dot(arr, iarr)

arr = np.array([[3, 2],
                [6, 4]])
print linalg.det(arr)
# print linalg.inv(arr)

print linalg.det(np.ones((3, 3)))

arr = np.arange(9).reshape((3, 3)) + np.diag([1, 0, 1])
print arr
uarr, spec, vharr = linalg.svd(arr)
print spec

sarr = np.diag(spec)
svd_mat = uarr.dot(sarr).dot(vharr)
print svd_mat
np.allclose(svd_mat, arr)
