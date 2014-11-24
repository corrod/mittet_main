from pylab import *
import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt('./out_diff/45pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]

#theta = omega *t
theta = 2.* np.pi* 2000. *x1
theta2 = theta -10

#radian = do * pi/ 180
radian = theta * np.pi/180

ax = plt.axes(polar=True)

ax.plot(theta,y1,'o',ls='-',ms=3,markevery=6)
# ax.plot(theta,y1,ls='-',markevery=6)

ylim(ymax=1e10)

plt.xlabel('45diff23every')

plt.savefig('polar45_23every.png')
# plt.show()
