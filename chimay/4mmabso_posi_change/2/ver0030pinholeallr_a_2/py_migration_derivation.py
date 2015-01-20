#coding : utf-8
import numpy as np
from pylab import *

left = 19
right = 113
cases = range(left,right)

data1 = np.loadtxt('./combine_pattern2_12diff.dat')
data2 = np.loadtxt('./combine_pattern2_23diff.dat')
data3 = np.loadtxt('./migration.dat')
data4 = np.loadtxt('./migration.dat')


#derivative of absolute
dy_stac = ( data1[:,0] - data1[:,1] ) / 2.
for i, case in enumerate(cases):
	y1 = data1[:,case-left+1]
	y2 = data1[:,case-left+2]
	dy = (y2-y1)/2.
	dy_stac = np.c_[dy_stac,dy]
np.savetxt('combine_pattern2_12derivation.dat',dy_stac)

a = len(dy_stac)
b = len(np.transpose(dy_stac))


fig = figure()
ax1 = fig.add_subplot(111)

ax1.set_xlabel('probe position')
ax1.set_ylabel('time')
ax1.set_title('deriviation of absolute')

x = arange(1+20,b+1+20) # x = arange(19,104+1)
y = arange(1,a+1) # y = arange(1,4351+1)
X,Y = meshgrid(x,y)

contourf(X,Y,dy_stac)
colorbar()
savefig('combine_pattern2_12derivation.png')
show()



#derivative of differential
dy_stac = ( data2[:,0] - data2[:,1] ) / 2.
for i, case in enumerate(cases):
	y1 = data2[:,case-left+1]
	y2 = data2[:,case-left+2]
	dy = (y2-y1)/2.
	dy_stac = np.c_[dy_stac,dy]
np.savetxt('combine_pattern2_23derivation.dat',dy_stac)

a = len(dy_stac)
b = len(np.transpose(dy_stac))


fig = figure()
ax1 = fig.add_subplot(111)

ax1.set_xlabel('probe position')
ax1.set_ylabel('time')
ax1.set_title('deriviation of differential')

x = arange(1,b+1) # x = arange(19,104+1)
y = arange(1,a+1) # y = arange(1,4351+1)
X,Y = meshgrid(x,y)

contourf(X,Y,dy_stac)
colorbar()
savefig('combine_pattern2_23derivation.png')
show()





left = 1
right = 94
cases = range(left,right)



######################################################################
# derivative of absolute migration
######################################################################

dy_stac = ( data3[:,0] - data3[:,1] ) / 2.
for i, case in enumerate(cases):
	y1 = data3[:,case]
	y2 = data3[:,case+1]
	dy = (y2-y1)/2.
	dy_stac = np.c_[dy_stac,dy]
np.savetxt('migration_abso_derivation.dat',dy_stac)

a = len(dy_stac)
b = len(np.transpose(dy_stac))


fig = figure()
ax1 = fig.add_subplot(111)

ax1.set_xlabel('probe position')
ax1.set_ylabel('time')
ax1.set_title('deriviation of absolute migration')

x = arange(1+20,b+1+20) # x = arange(19,104+1)
y = arange(1,a+1) # y = arange(1,4351+1)
X,Y = meshgrid(x,y)

contourf(X,Y,dy_stac)
colorbar()
savefig('migration_abso_derivation.png')
show()






# ######################################################################
# # derivative of differeantial migration
# ######################################################################

dy_stac = ( data4[:,0] - data4[:,1] ) / 2.
for i, case in enumerate(cases):
	y1 = data4[:,case]
	y2 = data4[:,case+1]
	dy = (y2-y1)/2.
	dy_stac = np.c_[dy_stac,dy]
np.savetxt('migration_diff_derivation.dat',dy_stac)

a = len(dy_stac)
b = len(np.transpose(dy_stac))


fig = figure()
ax1 = fig.add_subplot(111)

ax1.set_xlabel('probe position')
ax1.set_ylabel('time')
ax1.set_title('deriviation of differential migration')

x = arange(1,b+1) # x = arange(19,104+1)
y = arange(1,a+1) # y = arange(1,4351+1)
X,Y = meshgrid(x,y)

contourf(X,Y,dy_stac)
colorbar()
savefig('migration_diff_derivation.png')
show()