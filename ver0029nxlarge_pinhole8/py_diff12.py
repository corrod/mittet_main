# coding:utf-8

from pylab import *
import numpy as np
# import scipy as sp
import matplotlib.pyplot as plt

data1 = np.loadtxt('./out_diff/45pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
data2 = np.loadtxt('../ncver0027/out_diff/45pattern2_12diff.d')
x2 = data2[:,0]
y2 = data2[:,1]
plt.plot(x2,y2)
plt.xlabel('time[s]')
plt.ylabel('Differential Probe [A/m]')
plt.savefig('45diff12.png')
# plt.show()
