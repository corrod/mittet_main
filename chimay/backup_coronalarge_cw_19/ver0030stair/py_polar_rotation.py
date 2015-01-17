# coding:utf-8

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pylab import *

#define a list of markevery cases to plot
cases = [5,10,15,20,25
        	]


#define the figure size and grid layout properties
figsize = (7,6) #7,6 means 700,600points
cols = 3 #retsu
gs = gridspec.GridSpec(len(cases) // cols + 1, cols)


data1 = np.loadtxt('./out_diff/45pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]

#theta = omega *t
theta = 2.* np.pi* 2000. *x1
theta2 = theta -10

#radian = do * pi/ 180
radian = theta * np.pi/180

#plot each markevery case for polar plots
fig = plt.figure(num=4, figsize=figsize)
axpolar = []

for i, case in enumerate(cases):
    row = (i // cols)
    col = i % cols
    axpolar.append(fig.add_subplot(gs[row, col], polar = True))
    axpolar[-1].set_title('markevery=%s' % str(case))
    axpolar[-1].plot(theta+case, y1, 'o', ls='-', ms=3,  markevery=case)
# fig.tight_layout()

# ylim(ymax=1e10)

plt.xlabel('45diff23rotation')
plt.savefig('polar45_23rotation.png')
# plt.show()

