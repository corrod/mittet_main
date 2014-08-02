#coding utf-8
import numpy as np
a = np.array([0, 1, 2, 3])
print a
print a.ndim
print a.shape
print len(a)

b = np.array([[0, 1, 2], [3, 4, 5]])
print b

c = np.array([[[1], [2]], [[3], [4]]])
print c
print c.shape

d = np.arange(10)
print d

e = np.arange(1, 9, 2)  # arange:tousa
print e

f = np.linspace(0, 1, 6)  # n bunkatsu
print f

g = np.ones((3, 3))  # all one matrix
print g

h = np.zeros((3, 3))  # all zero matrix
print h

i = np.eye(3)  # identity matrix
print i

j = np.diag(np.array([1, 100, 3., -3]))  # diagonal matrix
print j

k = np.random.rand(3)
print k

l = np.ones((4, 4))
l[3, 1] = 6
l[2, 3] = 2
print l

m = np.diag(np.array([2., 3, 4, 5, 6]))
print m
