from __future__ import division, print_function

import numpy as np

from matplotlib.pyplot import *
from matplotlib.collections import LineCollection
import matplotlib.cbook as cbook

from pylab import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


rps=arange(41,61,2)

f = figure()
subplots_adjust(hspace=0.001)

# for rp in enumerate(rps):

# ax1 = subplot(911)
# data1 = np.loadtxt('./out_diff/41pattern2_23diff.d')
# x1 = data1[:,0]
# y1 = data1[:,1]
# plt.plot(x1,y1)

ax1 = subplot(911)
data1 = np.loadtxt('./out_diff/43pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
ylim(ymax=4e10)
ax1.set_yticklabels(['43'])
ax1.set_xticklabels([])

ax1 = subplot(912)
data1 = np.loadtxt('./out_diff/45pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
ylim(ymax=4e10)
ax1.set_yticklabels(['45'])
ax1.set_xticklabels([])

ax1 = subplot(913)
data1 = np.loadtxt('./out_diff/47pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
ylim(ymax=4e10)
ax1.set_yticklabels(['47'])
ax1.set_xticklabels([])

ax1 = subplot(914)
data1 = np.loadtxt('./out_diff/49pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
ylim(ymax=4e10)
ax1.set_yticklabels(['49'])
ax1.set_xticklabels([])

ax1 = subplot(915)
data1 = np.loadtxt('./out_diff/51pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
ylim(ymax=4e10)
ax1.set_yticklabels(['51'])
ax1.set_xticklabels([])

ax1 = subplot(916)
data1 = np.loadtxt('./out_diff/53pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
ylim(ymax=4e10)
ax1.set_yticklabels(['53'])
ax1.set_xticklabels([])

ax1 = subplot(917)
data1 = np.loadtxt('./out_diff/55pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
ylim(ymax=4e10)
ax1.set_yticklabels(['55'])
ax1.set_xticklabels([])

ax1 = subplot(918)
data1 = np.loadtxt('./out_diff/57pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
ylim(ymax=4e10)
ax1.set_yticklabels(['57'])
ax1.set_xticklabels([])

ax1 = subplot(919)
data1 = np.loadtxt('./out_diff/59pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
ax1.set_yticklabels(['59'])
ylim(ymax=4e10)
# ax1 = subplot(911)
# data1 = np.loadtxt('./out_diff/61pattern2_23diff.d')
# x1 = data1[:,0]
# y1 = data1[:,1]
# plt.plot(x1,y1)








# ax1 = subplot(1111)
# data1 = np.loadtxt('./out_diff/41pattern2_23diff.d')
# x1 = data1[:,0]
# y1 = data1[:,1]
# plt.plot(x1,y1)

# ax1 = subplot(1112)
# data1 = np.loadtxt('./out_diff/43pattern2_23diff.d')
# x1 = data1[:,0]
# y1 = data1[:,1]
# plt.plot(x1,y1)
# ax1 = subplot(1113)
# data1 = np.loadtxt('./out_diff/45pattern2_23diff.d')
# x1 = data1[:,0]
# y1 = data1[:,1]
# plt.plot(x1,y1)
# ax1 = subplot(1114)
# data1 = np.loadtxt('./out_diff/47pattern2_23diff.d')
# x1 = data1[:,0]
# y1 = data1[:,1]
# plt.plot(x1,y1)
# ax1 = subplot(1115)
# data1 = np.loadtxt('./out_diff/49pattern2_23diff.d')
# x1 = data1[:,0]
# y1 = data1[:,1]
# plt.plot(x1,y1)
# ax1 = subplot(1116)
# data1 = np.loadtxt('./out_diff/51pattern2_23diff.d')
# x1 = data1[:,0]
# y1 = data1[:,1]
# plt.plot(x1,y1)
# ax1 = subplot(1117)
# data1 = np.loadtxt('./out_diff/53pattern2_23diff.d')
# x1 = data1[:,0]
# y1 = data1[:,1]
# plt.plot(x1,y1)
# ax1 = subplot(1118)
# data1 = np.loadtxt('./out_diff/55pattern2_23diff.d')
# x1 = data1[:,0]
# y1 = data1[:,1]
# plt.plot(x1,y1)
# ax1 = subplot(1119)
# data1 = np.loadtxt('./out_diff/57pattern2_23diff.d')
# x1 = data1[:,0]
# y1 = data1[:,1]
# plt.plot(x1,y1)
# ax1 = subplot(11110)
# data1 = np.loadtxt('./out_diff/59pattern2_23diff.d')
# x1 = data1[:,0]
# y1 = data1[:,1]
# plt.plot(x1,y1)
# ax1 = subplot(11111)
# data1 = np.loadtxt('./out_diff/61pattern2_23diff.d')
# x1 = data1[:,0]
# y1 = data1[:,1]
# plt.plot(x1,y1)

savefig('timediffseries.png')
# show()

