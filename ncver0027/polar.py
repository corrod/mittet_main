from pylab import *
import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt('./out_diff/45pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
#theta = omega *t
theta = 2.* 2000.* 3.14159265359 *x1

ax = plt.axes(polar=True)

ax.scatter(theta,y1)
plt.savefig('ncpolar23.png')
# plt.show()