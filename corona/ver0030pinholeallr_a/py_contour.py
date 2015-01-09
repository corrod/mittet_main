#coding : utf-8
import numpy as np
from pylab import *
import matplotlib.pyplot as plt

#absolute
data1 = np.loadtxt('./combine_pattern2_12diff.dat')

x = arange(19,104+1)
y = arange(1,4351+1)
X,Y = meshgrid(x,y)

fig =figure()
ax1=fig.add_subplot(111)

ax1.set_xlabel('probe position')
ax1.set_ylabel('time')
ax1.set_title('absolute contour')

# contour(X,Y,data1) #no paste
contourf(X,Y,data1) #paste
colorbar()
savefig('combine_pattern2_12contour.png')
show()


# differential
data1 = np.loadtxt('./combine_pattern2_23diff.dat')

x = arange(19,104+1)
y = arange(1,4351+1)
X,Y = meshgrid(x,y)

fig =figure()
ax1=fig.add_subplot(111)

ax1.set_xlabel('probe position')
ax1.set_ylabel('time')
ax1.set_title('differential contour')

# contour(X,Y,data1) #no paste
contourf(X,Y,data1) #paste
colorbar()
savefig('combine_pattern2_23contour.png')
show()



# residual absolute
data1 = np.loadtxt('./combine_residual_abso.dat')

x = arange(19,104+1)
y = arange(1,4351+1)
X,Y = meshgrid(x,y)

fig =figure()
ax1=fig.add_subplot(111)

ax1.set_xlabel('probe position')
ax1.set_ylabel('time')
ax1.set_title('absolute residual contour')

# contour(X,Y,data1) #no paste
contourf(X,Y,data1) #paste
colorbar()
savefig('combine_residual_absocontour.png')
show()




# residual differential
data1 = np.loadtxt('./combine_residual_diff.dat')

x = arange(19,104+1)
y = arange(1,4351+1)
X,Y = meshgrid(x,y)

fig =figure()
ax1=fig.add_subplot(111)

ax1.set_xlabel('probe position')
ax1.set_ylabel('time')
ax1.set_title('differential residual contour')

# contour(X,Y,data1) #no paste
contourf(X,Y,data1) #paste
colorbar()
savefig('combine_residual_diffcontour.png')
show()