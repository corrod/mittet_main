# coding:utf-8
from pylab import *
import numpy as np
import matplotlib.pyplot as plt

cases = range(19,103)

for i, case in enumerate(cases):
	data1 = np.loadtxt("./out_diff/%spattern2_12diff.d" % str(case))
	t1 = data1[:,0]
	hz1 = data1[:,1]
	maxnumber=max(xrange(len(hz1)),key=lambda i:hz1[i])
	minnumber=min(xrange(len(hz1)),key=lambda i:hz1[i])
	print case,minnumber+1,min(hz1)
