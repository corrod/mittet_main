from pylab import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

cases = range(41,62,2)

# for i, case in enumerate(cases):
data1 = np.loadtxt('./out_diff/41pattern2_23diff.d')
t1 = data1[:,0]
hz1 = data1[:,1]
maxnumber=max(xrange(len(hz1)),key=lambda i:hz1[i])
minnumber=min(xrange(len(hz1)),key=lambda i:hz1[i])
print maxnumber+1, minnumber+1

