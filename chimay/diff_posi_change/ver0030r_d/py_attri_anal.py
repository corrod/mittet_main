# coding:utf-8
import numpy as np
import scipy.signal
from pylab import *
from math import *
from cmath import *

# data11 = np.loadtxt('combine_pattern2_12hilbert_r.dat')
# data22 = np.loadtxt('combine_pattern2_12hilbert_i.dat')
data11 = np.loadtxt('combine_pattern2_23hilbert_r.dat')
data22 = np.loadtxt('combine_pattern2_23hilbert_i.dat')
# data11 = np.loadtxt('migration_abso_hilbert_r.dat')
# data22 = np.loadtxt('migration_abso_hilbert_i.dat')
# data11 = np.loadtxt('migration_diff_hilbert_r.dat')
# data22 = np.loadtxt('migration_diff_hilbert_i.dat')

data1 = np.transpose(data11)
data2 = np.transpose(data22)
# data3 = np.transpose(data33)
# data4 = np.transpose(data44)
# data5 = np.transpose(data55)
# data6 = np.transpose(data66)


time_s = 8
time_e = 4350
# times = range(1,4351)

y_r = data1[time_s,:]
y_i = data2[time_s,:]
y_r_p = data1[time_s + 1, :]
y_r_m = data1[time_s - 1, :]
y_i_p = data2[time_s + 1, :]
y_i_m = data2[time_s - 1, :]
omega =( y_r * (y_i_p - y_i_m)/2. - y_i + (y_r_p - y_r_m)/2. ) / ( y_r**2 + y_i**2 )
omega_stac = omega

times = range(time_s + 1, time_e)

for i, time in enumerate(times):
	y_r = data1[time,:]
	y_i = data2[time,:]
	y_r_p = data1[time + 1, :]
	y_r_m = data1[time - 1, :]
	y_i_p = data2[time + 1, :]
	y_i_m = data2[time - 1, :]
	omega =( y_r * (y_i_p - y_i_m)/2. - y_i + (y_r_p - y_r_m)/2. ) / ( y_r**2 + y_i**2 )
	omega_stac = np.c_[omega_stac, omega]

omega_stac_t = np.transpose(omega_stac)
# print omega_stac

np.savetxt('instantaneous_freq_2_23.dat', omega_stac_t)

c = omega_stac_t
a = len(c)
b = len(np.transpose(c))
xlabel('Probe Position',fontsize=20)
ylabel('Time',fontsize=20)
title('Attribute Analysis instantaneous frequency differential')
x = arange(1,b+1)
y = arange(1,a+1)
X,Y = meshgrid(x,y)
contourf(X,Y,c)
colorbar()
savefig('instantaneous_freq_2_23.png')
plt.show()


plt.plot(omega_stac_t[:,50])
plt.title('instantaneous frequency one place')
savefig('instantaneous_freq_2_23_oneplace')
plt.show()

plt.plot(omega_stac_t[30,:])
plt.title('instantaneous frequency one time')
savefig('instantaneous_freq_2_23_onetime')
plt.show()













# data11 = np.loadtxt('combine_pattern2_12hilbert_r.dat')
# data22 = np.loadtxt('combine_pattern2_12hilbert_i.dat')
# data11 = np.loadtxt('combine_pattern2_23hilbert_r.dat')
# data22 = np.loadtxt('combine_pattern2_23hilbert_i.dat')
# data11 = np.loadtxt('migration_abso_hilbert_r.dat')
# data22 = np.loadtxt('migration_abso_hilbert_i.dat')
data11 = np.loadtxt('migration_diff_hilbert_r.dat')
data22 = np.loadtxt('migration_diff_hilbert_i.dat')

data1 = np.transpose(data11)
data2 = np.transpose(data22)
# data3 = np.transpose(data33)
# data4 = np.transpose(data44)
# data5 = np.transpose(data55)
# data6 = np.transpose(data66)


time_s = 8
time_e = 4350
# time_s = 2000
# time_e = 4350
# times = range(1,4351)

y_r = data1[time_s,:]
y_i = data2[time_s,:]
y_r_p = data1[time_s + 1, :]
y_r_m = data1[time_s - 1, :]
y_i_p = data2[time_s + 1, :]
y_i_m = data2[time_s - 1, :]
omega =( y_r * (y_i_p - y_i_m)/2. - y_i + (y_r_p - y_r_m)/2. ) / ( y_r**2 + y_i**2 )
omega_stac = omega

times = range(time_s + 1, time_e)

for i, time in enumerate(times):
	y_r = data1[time,:]
	y_i = data2[time,:]
	y_r_p = data1[time + 1, :]
	y_r_m = data1[time - 1, :]
	y_i_p = data2[time + 1, :]
	y_i_m = data2[time - 1, :]
	omega =( y_r * (y_i_p - y_i_m)/2. - y_i + (y_r_p - y_r_m)/2. ) / ( y_r**2 + y_i**2 )
	omega_stac = np.c_[omega_stac, omega]

omega_stac_t = np.transpose(omega_stac)
# print omega_stac

np.savetxt('instantaneous_freq_diffmigration.dat', omega_stac_t)

c = omega_stac_t
a = len(c)
b = len(np.transpose(c))
xlabel('Probe Position',fontsize=20)
ylabel('Time',fontsize=20)
title('Attribute Analysis instantaneous frequency differential migration')
x = arange(1,b+1)
y = arange(1,a+1)
X,Y = meshgrid(x,y)
xticks([22,42,62,82], ('40', '60', '80', '100'))
ylim(4351,0)
contourf(X,Y,c)
# colorbar()
savefig('instantaneous_freq_diffmigration.png')
plt.show()




# pick up ############################################
# c = omega_stac_t[3400:3600,38:58]
c = omega_stac_t[3400:3600,38:58]
a = len(c)
b = len(np.transpose(c))
xlabel('Probe Position',fontsize=20)
ylabel('Time',fontsize=20)
title('Attribute Analysis differential migration pickup')
x = arange(1,b+1)
y = arange(1,a+1)
X,Y = meshgrid(x,y)
xticks([1,5,10,15,20], ('55', '60', '65', '70', '75'))
# yticks([1,100,200,300,400], ('3300', '3400', '3500', '3600', '3700'))
yticks([1,100,200], ('3400', '3500', '3600'))
ylim(200,1)
contourf(X,Y,c)
colorbar()
savefig('instantaneous_freq_diffmigration_pickup.png')
plt.show()



#pickup one place
plt.plot(omega_stac_t[:,50])
plt.title('instantaneous frequency one place')
savefig('instantaneous_freq_diffmigration_oneplace')
plt.show()

#pickup one time
plt.plot(omega_stac_t[30,:])
plt.title('instantaneous frequency one time')
savefig('instantaneous_freq_diffmigration_onetime')
plt.show()


