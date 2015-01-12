# coding:utf-8
import numpy as np
from pylab import *

# time = range(0,1.737945116734669e-003,1.258468585615218e-006)
# np.savetxt('simulation_time.dat',time)
tt = 10

# t = tt * 1.258468585615218e-006
t = 10
rp = range(19,105)

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
ax1.legend()

ax2.set_title('extract differential at %s' % str(t))
ax2.set_xlabel('probe position')
ax2.set_ylabel('differential')
ax2.plot(rp,two_data,label='%s differential' % str(t))
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

ax3.set_title('extract absolute residual at %s' % str(t))
# ax1.set_xlabel('probe position')
ax3.set_ylabel('absolute residual')
ax3.plot(rp,one_data,label='%s absolute' % str(t))
ax3.legend()

ax4.set_title('extract differential residual at %s' % str(t))
ax4.set_xlabel('probe position')
ax4.set_ylabel('differential residual')
ax4.plot(rp,two_data,label='%s differential' % str(t))
ax4.legend()

savefig('%scombine_residual.png' % str(t))
show()