# coding:utf-8

from pylab import *
import numpy as np
import matplotlib.pyplot as plt

# import csv
# with open('coodinate.csv', 'wb') as f:
#     writer = csv.writer(f)
#     writer.writerows(someiterable)

cases = range(19,114)
f = open('./pycordinate.dat','w')
# for case in enumerate(cases):
# 	print >> f, case
# f.close

# for ix in range(19,104):
	# print >> f, ix, ","
g=range(19,104)
print >> f, range(19,104)
f.close