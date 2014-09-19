from pylab import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


f = figure()
subplots_adjust(hspace=0.001)

ax1 = subplot(611)
data1 = np.loadtxt('./hz1000.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)

ax2 = subplot(612)
data2 = np.loadtxt('./hz1010.d')
x2 = data2[:,0]
y2 = data2[:,1]
plt.plot(x2,y2)

ax3 = subplot(613)
data3 = np.loadtxt('./hz1020.d')
x3 = data3[:,0]
y3 = data3[:,1]
plt.plot(x3,y3)

ax4 = subplot(614)
data4 = np.loadtxt('./hz1030.d')
x4 = data4[:,0]
y4 = data4[:,1]
plt.plot(x4,y4)

ax5 = subplot(615)
data5 = np.loadtxt('./hz1040.d')
x5 = data5[:,0]
y5 = data5[:,1]
plt.plot(x5,y5)

ax6 = subplot(616)
data6 = np.loadtxt('./hz1050.d')
x6 = data6[:,0]
y6 = data6[:,1]
plt.plot(x6,y6)

show()



# from pylab import *

# t = arange(0.0, 2.0, 0.01)

# s1 = sin(2*pi*t)
# s2 = exp(-t)
# s3 = s1*s2

# # axes rect in relative 0,1 coords left, bottom, width, height.  Turn
# # off xtick labels on all but the lower plot


# f = figure()
# subplots_adjust(hspace=0.001)


# ax1 = subplot(311)
# ax1.plot(t,s1)
# yticks(arange(-0.9, 1.0, 0.4))
# ylim(-1,1)

# ax2 = subplot(312, sharex=ax1)
# ax2.plot(t,s2)
# yticks(arange(0.1, 1.0,  0.2))
# ylim(0,1)

# ax3 = subplot(313, sharex=ax1)
# ax3.plot(t,s3)
# yticks(arange(-0.9, 1.0, 0.4))
# ylim(-1,1)

# xticklabels = ax1.get_xticklabels()+ax2.get_xticklabels()
# setp(xticklabels, visible=False)

# show()
