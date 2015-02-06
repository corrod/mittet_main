# coding:utf-8
import numpy as np
import scipy.signal
from pylab import *
from math import *
from cmath import *

data4 = np.loadtxt('./migration.dat')



#cross correlation
left = 19
right = 112
cases = range(left,right)

y4 = data4[:,0]
y42= data4[:,1]

y_corr4 = scipy.signal.correlate(y4,y42,"same")
y_corr_stac4 = y_corr4

for i, case in enumerate(cases):
	y4 = data4[:,case-left+1]
	y42 = data4[:,case-left+2]
	y_corr4 = scipy.signal.correlate(y4,y42,"same")
	y_corr_stac4  = np.c_[y_corr_stac4 ,y_corr4]

np.savetxt('migration_correlation.dat',y_corr_stac4)



#plot cross correlation
a = len(y_corr_stac4)
b = len(np.transpose(y_corr_stac4))

fig = figure()
ax1 = fig.add_subplot(121)
title('Cross Correlation of migration')
ax1.set_xlabel('Probe Position',fontsize=20)
ax1.set_ylabel('Time',fontsize=20)
# ax1.set_title('')

x = arange(1,b+1)
y = arange(1,a+1)
xticks([22,42,62,82], ('40', '60', '80', '100'))
# yticks([500,1000,1500,2000,2500,3000,3500,4000], (''))
X,Y = meshgrid(x,y)

contourf(X,Y,y_corr_stac4)
# colorbar()
# savefig('migration_corr.png')
# show()




#hilbert transform
cases = range(1,4351)

y3 = y_corr_stac4[0,:]
y_hilbert3 = scipy.signal.hilbert(y3)
y_hilbert_stac3 = y_hilbert3
envelope_stac3 = (real(y_hilbert3)**2 + imag(y_hilbert3)**2)**0.5
for i, case in enumerate(cases):
	y3 = y_corr_stac4[case,:]
	y_hilbert3 = scipy.signal.hilbert(y3)
	y_hilbert_stac3 = np.c_[y_hilbert_stac3,y_hilbert3]
	envelope3 = (real(y_hilbert3)**2 + imag(y_hilbert3)**2)**0.5
	envelope_stac3 = np.c_[envelope_stac3,envelope3]

y_hilbert_stac3_r = real(y_hilbert_stac3)
y_hilbert_stac3_i = imag(y_hilbert_stac3)

np.savetxt('migration_corr_hilbert_r.dat',y_hilbert_stac3_r)
np.savetxt('migration_corr_hilbert_i.dat',y_hilbert_stac3_i)
np.savetxt('migration_corr_hilbert_envelope',envelope_stac3)


# #hilbert transform
# cases = range(1,94)
# # cases = range(1,4351)

# y3 = y_corr_stac4[:,0]
# y_hilbert3 = scipy.signal.hilbert(y3)
# y_hilbert_stac3 = y_hilbert3
# envelope_stac3 = (real(y_hilbert3)**2 + imag(y_hilbert3)**2)**0.5
# for i, case in enumerate(cases):
# 	y3 = y_corr_stac4[:,case]
# 	y_hilbert3 = scipy.signal.hilbert(y3)
# 	y_hilbert_stac3 = np.c_[y_hilbert_stac3,y_hilbert3]
# 	envelope3 = (real(y_hilbert3)**2 + imag(y_hilbert3)**2)**0.5
# 	envelope_stac3 = np.c_[envelope_stac3,envelope3]

# y_hilbert_stac3_r = real(y_hilbert_stac3)
# y_hilbert_stac3_i = imag(y_hilbert_stac3)

# np.savetxt('migration_corr_hilbert_r.dat',y_hilbert_stac3_r)
# np.savetxt('migration_corr_hilbert_i.dat',y_hilbert_stac3_i)
# np.savetxt('migration_corr_hilbert_envelope',envelope_stac3)





#plot envelope
c = np.transpose(envelope_stac3)
a = len(c)
b = len(np.transpose(c))

ax1=subplot(122)
# ax1=subplot(121)
xlabel('Probe Position',fontsize=20)
ylabel('Time',fontsize=20)
title('Envelope of Cross Correlation')
x = arange(1,b+1)
y = arange(1,a+1)
xticks([22,42,62,82], ('40', '60', '80', '100'))

ax1.set_yticklabels([])
X,Y = meshgrid(x,y)
contourf(X,Y,c)
# colorbar()
savefig('migration_corr_envelope.png')
show()


#####################################################
# pickup max time minimum time
# tabun yarikata machigatteru
#####################################################
cases = range(0,94)
time_max=0
time_min=0
# ev = np.transpose(envelope_stac3)
#plot max time
for i, case in enumerate(cases):
	hz1 = envelope_stac3[:,case]
	# hz1 = ev[:,case]
	maxnumber=max(xrange(len(hz1)),key=lambda i:hz1[i])
	minnumber=min(xrange(len(hz1)),key=lambda i:hz1[i])
	time_max = np.c_[time_max,maxnumber]
	time_min = np.c_[time_min,minnumber]

	# print case,maxnumber+1,max(hz1)
# print time_max, time_min
time_max_t = np.transpose(time_max)
time_min_t = np.transpose(time_min)


plot(time_max_t[1:])
# title( "max time")
xlabel('Probe Position',fontsize=20)
ylabel('Time',fontsize=20)
xticks([23,43,63,83], ('40', '60', '80', '100'))
savefig('migration_corr_envelope_maxtimefalse.png')
show()

# plot(time_min_t[1:])
plot(time_min_t[1:91])
# title( "min time")
xlabel('Probe Position',fontsize=20)
ylabel('Time',fontsize=20)
xticks([23,43,63,83], ('40', '60', '80', '100'))
savefig('migration_corr_envelope_mintimefalse.png')
show()


###########################################################
# # pickup max time minimum time
#   tabun kotti ga hontou
########################################################
cases = range(0,94)
time_max=0
time_min=0
ev = np.transpose(envelope_stac3)
#plot max time
for i, case in enumerate(cases):
	# hz1 = envelope_stac3[case,:]
	hz1 = ev[:,case]
	maxnumber=max(xrange(len(hz1)),key=lambda i:hz1[i])
	minnumber=min(xrange(len(hz1)),key=lambda i:hz1[i])
	time_max = np.c_[time_max,maxnumber]
	time_min = np.c_[time_min,minnumber]

	# print case,maxnumber+1,max(hz1)
# print time_max, time_min
time_max_t = np.transpose(time_max)
time_min_t = np.transpose(time_min)


plot(time_max_t[1:])
# title( "max time")
xlabel('Probe Position',fontsize=20)
ylabel('Time',fontsize=20)
xticks([23,43,63,83], ('40', '60', '80', '100'))
savefig('migration_corr_envelope_maxtimetrue.png')
show()

plot(time_min_t[1:])
# title( "min time")
xlabel('Probe Position',fontsize=20)
ylabel('Time',fontsize=20)
xticks([23,43,63,83], ('40', '60', '80', '100'))
savefig('migration_corr_envelope_mintimetrue.png')
show()



# how to cross correlate

# x = [1.0, 2.0, 3.0, 4.]
# y = [1.0, 2.0, 3.0, 4.]
# y = [120.0, 130.0, 90.0]
# x=[1, 2, 3]
# y=[0, 1, 0.5]
# cor_varid = scipy.signal.correlate(x,y,mode='valid')
# cor_same = scipy.signal.correlate(x,y,mode='same')
# cor_full = scipy.signal.correlate(x,y,mode='full')
# print cor_varid
# print cor_same
# print cor_full

# cor_varid_np = np.correlate(x,y,mode='valid')
# cor_same_np = np.correlate(x,y,mode='same')
# cor_full_np = np.correlate(x,y,mode='full')
# print cor_varid_np
# print cor_same_np
# print cor_full_np

