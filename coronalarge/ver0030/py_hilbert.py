import numpy as np
import scipy as sp
from pylab import *

# left = 19
# right = 113
# cases = range(left,right)

data1 = np.loadtxt('./combine_pattern2_12diff.dat')
data2 = np.loadtxt('./combine_pattern2_23diff.dat')
data3 = np.loadtxt('./combine_residual_abso.dat')
data4 = np.loadtxt('./combine_residual_diff.dat')

cases = range(1,4351)
#hilbert transform of absolute
for i, case in enumerato(cases):
	y[case,:] = data1[case,:]
	y_hilbert = sp.signal.hilbert(y)
np.savetxt('combine_pattern2_12hilbert.dat',y_hilbert)

#hilbert transform of differential
#hilbert transform of absolute residual
#hilbert transform of differential residual