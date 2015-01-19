# coding:utf-8
import numpy as np
import scipy.signal
from pylab import *
from math import *
from cmath import *


data1 = np.loadtxt('./combine_pattern2_12diff.dat')
data2 = np.loadtxt('./combine_pattern2_23diff.dat')
data3 = np.loadtxt('./migration.dat')
data4 = np.loadtxt('./migration.dat')

cases = range(1,4351)




######################################################################
#hilbert transform of absolute
######################################################################
y1 = data1[0,:]
y_hilbert1 = scipy.signal.hilbert(y1)
y_hilbert_stac1 = y_hilbert1
envelope_stac1 = (real(y_hilbert1)**2 + imag(y_hilbert1)**2)**0.5
for i, case in enumerate(cases):
	y1 = data1[case,:]
	y_hilbert1 = scipy.signal.hilbert(y1)
	y_hilbert_stac1 = np.c_[y_hilbert_stac1,y_hilbert1]
	envelope1 = (real(y_hilbert1)**2 + imag(y_hilbert1)**2)**0.5
	envelope_stac1 = np.c_[envelope_stac1,envelope1]
y_hilbert_stac1_r = real(y_hilbert_stac1)
y_hilbert_stac1_i = imag(y_hilbert_stac1)

np.savetxt('combine_pattern2_12hilbert_r.dat',y_hilbert_stac1_r)
np.savetxt('combine_pattern2_12hilbert_i.dat',y_hilbert_stac1_i)
np.savetxt('combine_pattern2_12hilbert_envelope',envelope_stac1)

fig = figure(figsize=(15,5))


c = np.transpose(y_hilbert_stac1_r)
a = len(c)
b = len(np.transpose(c))

subplot(131)
xlabel('probe position')
ylabel('time')
title('hilbert real of absolute')
x = arange(1+19,b+1+19)
y = arange(1,a+1)
X,Y = meshgrid(x,y)
contourf(X,Y,c)
colorbar()


c = np.transpose(y_hilbert_stac1_i)
a = len(c)
b = len(np.transpose(c))

subplot(132)
xlabel('probe position')
ylabel('time')
title('hilbert imag of absolute')
x = arange(1+19,b+1+19)
y = arange(1,a+1)
X,Y = meshgrid(x,y)
contourf(X,Y,c)
colorbar()
# savefig('combine_pattern_2_12_hilbert.png')
# show()

c = np.transpose(envelope_stac1)
a = len(c)
b = len(np.transpose(c))

subplot(133)
xlabel('probe position')
ylabel('time')
title('hilbert envelope of absolute')
x = arange(1+19,b+1+19)
y = arange(1,a+1)
X,Y = meshgrid(x,y)
contourf(X,Y,c)
colorbar()
savefig('combine_pattern_2_12_hilbert_envelope.png')
show()





######################################################################
#hilbert transform of differential
######################################################################

y2 = data2[0,:]
y_hilbert2 = scipy.signal.hilbert(y2)
y_hilbert_stac2 = y_hilbert2
envelope_stac2 = (real(y_hilbert2)**2 + imag(y_hilbert2)**2)**0.5
for i, case in enumerate(cases):
	y2 = data2[case,:]
	y_hilbert2 = scipy.signal.hilbert(y2)
	y_hilbert_stac2 = np.c_[y_hilbert_stac2,y_hilbert2]
	envelope2 = (real(y_hilbert2)**2 + imag(y_hilbert2)**2)**0.5
	envelope_stac2 = np.c_[envelope_stac2,envelope2]
y_hilbert_stac2_r = real(y_hilbert_stac2)
y_hilbert_stac2_i = imag(y_hilbert_stac2)
np.savetxt('combine_pattern2_23hilbert_r.dat',y_hilbert_stac2_r)
np.savetxt('combine_pattern2_23hilbert_i.dat',y_hilbert_stac2_i)
np.savetxt('combine_pattern2_23hilbert_envelope',envelope_stac2)


fig = figure(figsize=(15,5))


c = np.transpose(y_hilbert_stac2_r)
a = len(c)
b = len(np.transpose(c))

subplot(131)
xlabel('probe position')
ylabel('time')
title('hilbert real of differential')
x = arange(1,b+1)
y = arange(1,a+1)
X,Y = meshgrid(x,y)
contourf(X,Y,c)
colorbar()


c = np.transpose(y_hilbert_stac2_i)
a = len(c)
b = len(np.transpose(c))

subplot(132)
xlabel('probe position')
ylabel('time')
title('hilbert imag of differential')
x = arange(1,b+1)
y = arange(1,a+1)
X,Y = meshgrid(x,y)
contourf(X,Y,c)
colorbar()
# savefig('combine_pattern_2_23_hilbert.png')
# show()



c = np.transpose(envelope_stac2)
a = len(c)
b = len(np.transpose(c))

subplot(133)
xlabel('probe position')
ylabel('time')
title('hilbert envelope of differential')
x = arange(1,b+1)
y = arange(1,a+1)
X,Y = meshgrid(x,y)
contourf(X,Y,c)
colorbar()
savefig('combine_pattern_2_23_hilbert_envelope.png')
show()






######################################################################
#hilbert transform of absolute migration
######################################################################

y3 = data3[0,:]
y_hilbert3 = scipy.signal.hilbert(y3)
y_hilbert_stac3 = y_hilbert3
envelope_stac3 = (real(y_hilbert3)**2 + imag(y_hilbert3)**2)**0.5
for i, case in enumerate(cases):
	y3 = data3[case,:]
	y_hilbert3 = scipy.signal.hilbert(y3)
	y_hilbert_stac3 = np.c_[y_hilbert_stac3,y_hilbert3]
	envelope3 = (real(y_hilbert3)**2 + imag(y_hilbert3)**2)**0.5
	envelope_stac3 = np.c_[envelope_stac3,envelope3]
y_hilbert_stac3_r = real(y_hilbert_stac3)
y_hilbert_stac3_i = imag(y_hilbert_stac3)
np.savetxt('migration_abso_hilbert_r.dat',y_hilbert_stac3_r)
np.savetxt('migration_abso_hilbert_i.dat',y_hilbert_stac3_i)
np.savetxt('migration_abso_hilbert_envelope',envelope_stac3)


fig = figure(figsize=(15,5))


c = np.transpose(y_hilbert_stac3_r)
a = len(c)
b = len(np.transpose(c))

subplot(131)
xlabel('probe position')
ylabel('time')
title('hilbert real of absolute migration')
x = arange(1,b+1)
y = arange(1,a+1)
X,Y = meshgrid(x,y)
contourf(X,Y,c)
colorbar()


c = np.transpose(y_hilbert_stac3_i)
a = len(c)
b = len(np.transpose(c))

subplot(132)
xlabel('probe position')
ylabel('time')
title('hilbert imag of absolute migration')
x = arange(1,b+1)
y = arange(1,a+1)
X,Y = meshgrid(x,y)
contourf(X,Y,c)
colorbar()
# savefig('migration_abso_hilbert.png')
# show()



c = np.transpose(envelope_stac3)
a = len(c)
b = len(np.transpose(c))

subplot(133)
xlabel('probe position')
ylabel('time')
title('hilbert envelope of absolute migration')
x = arange(1,b+1)
y = arange(1,a+1)
X,Y = meshgrid(x,y)
contourf(X,Y,c)
colorbar()
savefig('migration_abso_hilbert_envelope.png')
show()


######################################################################
#hilbert transform of differential migration
######################################################################

y4 = data4[0,:]
y_hilbert4 = scipy.signal.hilbert(y4)
y_hilbert_stac4 = y_hilbert4
envelope_stac4 = (real(y_hilbert4)**2 + imag(y_hilbert4)**2)**0.5
for i, case in enumerate(cases):
	y4 = data4[case,:]
	y_hilbert4 = scipy.signal.hilbert(y4)
	y_hilbert_stac4 = np.c_[y_hilbert_stac4,y_hilbert4]
	envelope4 = (real(y_hilbert4)**2 + imag(y_hilbert4)**2)**0.5
	envelope_stac4 = np.c_[envelope_stac4,envelope4]
y_hilbert_stac4_r = real(y_hilbert_stac4)
y_hilbert_stac4_i = imag(y_hilbert_stac4)
np.savetxt('migration_diff_hilbert_r.dat',y_hilbert_stac4_r)
np.savetxt('migration_diff_hilbert_i.dat',y_hilbert_stac4_i)
np.savetxt('migration_diff_hilbert_envelope',envelope_stac4)


fig = figure(figsize=(15,5))


c = np.transpose(y_hilbert_stac4_r)
a = len(c)
b = len(np.transpose(c))

subplot(131)
xlabel('probe position')
ylabel('time')
title('hilbert real of differential migration')
x = arange(1,b+1)
y = arange(1,a+1)
X,Y = meshgrid(x,y)
contourf(X,Y,c)
colorbar()


c = np.transpose(y_hilbert_stac4_i)
a = len(c)
b = len(np.transpose(c))

subplot(132)
xlabel('probe position')
ylabel('time')
title('hilbert imag of differential migration')
x = arange(1,b+1)
y = arange(1,a+1)
X,Y = meshgrid(x,y)
contourf(X,Y,c)
colorbar()
# savefig('migration_diff_hilbert.png')
# show()



c = np.transpose(envelope_stac4)
a = len(c)
b = len(np.transpose(c))

subplot(133)
xlabel('probe position')
ylabel('time')
title('hilbert envelope of differential migration')
x = arange(1,b+1)
y = arange(1,a+1)
X,Y = meshgrid(x,y)
contourf(X,Y,c)
colorbar()
savefig('migration_diff_hilbert_envelope.png')
show()