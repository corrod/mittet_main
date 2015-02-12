# coding:utf-8
import numpy as np
from pylab import *

# time = range(0,1.737945116734669e-003,1.258468585615218e-006)
# np.savetxt('simulation_time.dat',time)

t = 10

rp = range(19,115)

fig = figure()

ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

data1 = np.loadtxt('./combine_pattern2_12diff.dat')
data2 = np.loadtxt('./combine_pattern2_23diff.dat')
# data1 = np.loadtxt('./combine_residual_abso.dat')
# data2 = np.loadtxt('./combine_residual_diff.dat')



one_data = data1[t,:]
two_data = data2[t,:]

ax1.set_title('extract absolute at %s' % str(t))
# ax1.set_xlabel('probe position')
ax1.set_ylabel('absolute')
ax1.plot(rp,one_data,label='%s absolute' % str(t))
ax1.grid(True)
ax1.legend()

ax2.set_title('extract differential at %s' % str(t))
ax2.set_xlabel('probe position')
ax2.set_ylabel('differential')
ax2.plot(rp,two_data,label='%s differential' % str(t))
ax2.grid(True)
ax2.legend()

savefig('%scombine_pattern.png' % str(t))
show()


fig=figure()


ax3 = fig.add_subplot(211)
ax4 = fig.add_subplot(212)

# data1 = np.loadtxt('./combine_pattern2_12diff.dat')
# data2 = np.loadtxt('./combine_pattern2_23diff.dat')
data1 = np.loadtxt('./combine_residual_abso.dat')
data2 = np.loadtxt('./combine_residual_diff.dat')



one_data = data1[t,:]
two_data = data2[t,:]

# ax3.set_title('extract absolute residual at %s' % str(t))
# # ax1.set_xlabel('probe position')
# ax3.set_ylabel('absolute residual')
# ax3.plot(rp,one_data,label='%s absolute' % str(t))
# ax3.grid(True)
# ax3.legend()

ax4.set_title('extract differential residual at one time')
ax4.set_xlabel('Probe Position',fontsize=20)
ax4.set_ylabel('Residual [A/m]',fontsize=20)
ax4.plot([20,20],[-4e-7,4e-7], color='red', linewidth=1, linestyle="--")
ax4.plot([30,30],[-4e-7,4e-7], color='red', linewidth=1, linestyle="--")
ax4.plot([39,39],[-4e-7,4e-7], color='red', linewidth=1, linestyle="--")
ax4.plot([41,41],[-4e-7,4e-7], color='red', linewidth=1, linestyle="--")
ax4.plot([53,53],[-4e-7,4e-7], color='red', linewidth=1, linestyle="--")
ax4.plot([57,57],[-4e-7,4e-7], color='red', linewidth=1, linestyle="--")
ax4.plot([67,67],[-4e-7,4e-7], color='red', linewidth=1, linestyle="--")
ax4.plot([73,73],[-4e-7,4e-7], color='red', linewidth=1, linestyle="--")
ax4.plot([93,93],[-4e-7,4e-7], color='red', linewidth=1, linestyle="--")
ax4.plot([107,107],[-4e-7,4e-7], color='red', linewidth=1, linestyle="--")
ax4.plot(rp,two_data,color='blue',label='diff')
ax4.grid(True)
ax4.legend()

savefig('%scombine_residual.png' % str(t))
show()