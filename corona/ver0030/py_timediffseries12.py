# coding:utf-8
from __future__ import division, print_function

import numpy as np

from matplotlib.pyplot import *
from matplotlib.collections import LineCollection
import matplotlib.cbook as cbook

from pylab import *
import numpy as np
# import scipy as sp
import matplotlib.pyplot as plt


rps=arange(41,61,2)

f = figure()
subplots_adjust(hspace=0.001)


# for rp in enumerate(rps):

# ax1 = subplot(911)
# data1 = np.loadtxt('./out_diff/41pattern2_12diff.d')
# x1 = data1[:,0]
# y1 = data1[:,1]
# plt.plot(x1,y1)
ax1 = subplot(911)
ax1.set_title('absolute method')
data1 = np.loadtxt('./out_diff/43pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
# ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('43')

ax1 = subplot(912)
data1 = np.loadtxt('./out_diff/45pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('45')

ax1 = subplot(913)
data1 = np.loadtxt('./out_diff/47pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('47')

ax1 = subplot(914)
data1 = np.loadtxt('./out_diff/49pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('49')

ax1 = subplot(915)
data1 = np.loadtxt('./out_diff/51pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('51')

ax1 = subplot(916)
data1 = np.loadtxt('./out_diff/53pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('53')

ax1 = subplot(917)
data1 = np.loadtxt('./out_diff/55pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('55')

ax1 = subplot(918)
data1 = np.loadtxt('./out_diff/57pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('57')

ax1 = subplot(919)
data1 = np.loadtxt('./out_diff/59pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
ax1.set_yticklabels([])
# ylim(ymax=2e11)
ax1.set_ylabel('59')

# ax1 = subplot(911)
# data1 = np.loadtxt('./out_diff/61pattern2_12diff.d')
# x1 = data1[:,0]
# y1 = data1[:,1]
# plt.plot(x1,y1)

plt.xlabel('time [s]')



savefig('timediffseries12.png')
show()

