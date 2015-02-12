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
# ax1.set_xlabel('Probe Position')
ax1.set_ylabel('Absolute [A/m]')
ax1.plot(rp,one_data,label='%s absolute' % str(t))
ax1.grid(True)
ax1.legend()

ax2.set_title('extract differential at %s' % str(t))
ax2.set_xlabel('Probe Position')
ax2.set_ylabel('Differential [A/m]')
ax2.plot(rp,two_data,label='%s differential' % str(t))
ax2.grid(True)
ax2.legend()

savefig('%scombine_pattern.png' % str(t))
show()


fig=figure()
# fig = figure(figsize=(10,5))


ax3 = subplot(211)
# ax3 = subplot(111)

# ax1 = fig.add_subplot(111)


# data1 = np.loadtxt('./combine_pattern2_12diff.dat')
# data2 = np.loadtxt('./combine_pattern2_23diff.dat')
data1 = np.loadtxt('./combine_residual_abso.dat')
data2 = np.loadtxt('./combine_residual_diff.dat')



one_data = data1[t,:]
two_data = data2[t,:]

ax3.set_title('extract absolute residual at one time')
ax3.set_xlabel('Probe Position',fontsize=20)
ax3.set_ylabel('Residual [A/m]',fontsize=20)
ax3.plot([64,64],[0,4e-7],color='red' , linewidth=1, linestyle="--")
# ax3.plot([65,65],[0,4e-7], color='red', linewidth=3, linestyle="--")
# ax3.plot([66,66],[0,4e-7], color='red', linewidth=3, linestyle="--")
# ax3.plot([67,67],[0,4e-7], color='red', linewidth=3, linestyle="--")
ax3.plot([68,68],[0,4e-7], color='red', linewidth=1, linestyle="--")
ax3.plot(rp,one_data,color='blue',label='%s absolute' % str(t))
ax3.grid(True)
ax3.legend()
# savefig('%scombine_absoresidual.png' % str(t))
# show()

# # fig = figure(figsize=(9,4))
# ax4 = subplot(212)
# # ax4 = subplot(111)
# ax4.set_title('extract differential residual at one time')
# ax4.set_xlabel('Probe Position',fontsize=20)
# ax4.set_ylabel('Residual [A/m]',fontsize=20)
# ax4.plot(rp,two_data,label='%s differential' % str(t))
# ax4.grid(True)
# ax4.legend()

savefig('%scombine_residual.png' % str(t))
show()








# # coding:utf-8
# import numpy as np
# from pylab import *

# # time = range(0,1.737945116734669e-003,1.258468585615218e-006)
# # np.savetxt('simulation_time.dat',time)

# t = 10

# rp = range(19,115)

# fig = figure()

# ax1 = fig.add_subplot(211)
# ax2 = fig.add_subplot(212)

# data1 = np.loadtxt('./combine_pattern2_12diff.dat')
# data2 = np.loadtxt('./combine_pattern2_23diff.dat')
# # data1 = np.loadtxt('./combine_residual_abso.dat')
# # data2 = np.loadtxt('./combine_residual_diff.dat')



# one_data = data1[t,:]
# two_data = data2[t,:]

# ax1.set_title('extract absolute at %s' % str(t))
# # ax1.set_xlabel('Probe Position')
# ax1.set_ylabel('Absolute [A/m]')
# ax1.plot(rp,one_data,label='%s absolute' % str(t))
# ax1.grid(True)
# ax1.legend()

# ax2.set_title('extract differential at %s' % str(t))
# ax2.set_xlabel('Probe Position')
# ax2.set_ylabel('Differential [A/m]')
# ax2.plot(rp,two_data,label='%s differential' % str(t))
# ax2.grid(True)
# ax2.legend()

# savefig('%scombine_pattern.png' % str(t))
# show()


# # fig=figure()
# fig = figure()


# ax3 = subplot(111)

# # ax1 = fig.add_subplot(111)


# # data1 = np.loadtxt('./combine_pattern2_12diff.dat')
# # data2 = np.loadtxt('./combine_pattern2_23diff.dat')
# data1 = np.loadtxt('./combine_residual_abso.dat')
# data2 = np.loadtxt('./combine_residual_diff.dat')



# one_data = data1[t,:]
# two_data = data2[t,:]

# ax3.set_title('extract absolute residual at one time')
# ax3.set_xlabel('Probe Position',fontsize=20)
# ax3.set_ylabel('Residual [A/m]',fontsize=20)
# ax3.plot(rp,one_data,label='%s absolute' % str(t))
# ax3.grid(True)
# ax3.legend()
# savefig('%scombine_absoresidual.png' % str(t))
# show()

# fig = figure()
# ax4 = subplot(111)
# ax4.set_title('extract differential residual at one time')
# ax4.set_xlabel('Probe Position')
# ax4.set_ylabel('Residual [A/m]')
# ax4.plot(rp,two_data,label='%s differential' % str(t))
# ax4.grid(True)
# ax4.legend()

# savefig('%scombine_residual.png' % str(t))
# show()