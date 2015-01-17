# coding:utf-8
from matplotlib.pyplot import *
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
# subplots_adjust(hspace=0.001)

fig=figure()

ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)


data1 = np.loadtxt('./abso_max.d')
posi1 = data1[:,0]
time1 = data1[:,1]
max1  = data1[:,2]

ax1.set_title('extract max at each place')
# ax1.set_xlabel('shot position (left)')
# ax1.set_xticks([])
ax1.set_ylabel('absolute max')
ax1.plot(posi1,max1,label='max')
ax1.legend()


data2 = np.loadtxt('./abso_min.d')
posi2 = data2[:,0]
time2 = data2[:,1]
min2  = data2[:,2]

ax2.set_title('extract min at each place')
ax2.set_xlabel('shot position x (2)')
ax2.set_ylabel('absolute min')
ax2.plot(posi2,min2,label='min')
ax2.legend()

# fig.tight_layout()
savefig('absomaxminplot.png')
show() # after savefig