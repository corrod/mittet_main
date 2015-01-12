# coding:utf-8
import numpy as np
import scipy.signal
from pylab import *


data1 = np.loadtxt('./combine_pattern2_12diff.dat')
data2 = np.loadtxt('./combine_pattern2_23diff.dat')
data3 = np.loadtxt('./combine_residual_abso.dat')
data4 = np.loadtxt('./combine_residual_diff.dat')

cases = range(1,4351)

#hilbert transform of absolute
y1 = data1[0,:]
y_hilbert1 = scipy.signal.hilbert(y1)
y_hilbert_stac1 = y_hilbert1
for i, case in enumerate(cases):
	y1 = data1[case,:]
	y_hilbert1 = scipy.signal.hilbert(y1)
	y_hilbert_stac1 = np.c_[y_hilbert_stac1,y_hilbert1]
np.savetxt('combine_pattern2_12hilbert.dat',y_hilbert_stac1)


#hilbert transform of differential
y2 = data2[0,:]
y_hilbert2 = scipy.signal.hilbert(y2)
y_hilbert_stac2 = y_hilbert2
for i, case in enumerate(cases):
	y2 = data2[case,:]
	y_hilbert2 = scipy.signal.hilbert(y2)
	y_hilbert_stac2 = np.c_[y_hilbert_stac2,y_hilbert2]
np.savetxt('combine_pattern2_23hilbert.dat',y_hilbert_stac2)


#hilbert transform of absolute residual
y3 = data3[0,:]
y_hilbert3 = scipy.signal.hilbert(y3)
y_hilbert_stac3 = y_hilbert3
for i, case in enumerate(cases):
	y3 = data3[case,:]
	y_hilbert3 = scipy.signal.hilbert(y3)
	y_hilbert_stac3 = np.c_[y_hilbert_stac3,y_hilbert3]
np.savetxt('combine_residual_abso_hilbert.dat',y_hilbert_stac3)


#hilbert transform of differential residual
y4 = data4[0,:]
y_hilbert4 = scipy.signal.hilbert(y4)
y_hilbert_stac4 = y_hilbert4
for i, case in enumerate(cases):
	y4 = data4[case,:]
	y_hilbert4 = scipy.signal.hilbert(y4)
	y_hilbert_stac4 = np.c_[y_hilbert_stac4,y_hilbert4]
np.savetxt('combine_residual_diff_hilbert.dat',y_hilbert_stac4)
