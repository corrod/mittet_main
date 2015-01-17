# coding:utf-8

import numpy as np
from pylab import *

inputdata = np.loadtxt('./signal.d')
t_in = inputdata[:,0]
hz_r_in = inputdata[:,1]
hz_i_in = inputdata[:,2]

cases = range(19,103)

for i, case in enumerate(cases):
	data1 = np.loadtxt("./out_diff/%spattern2_12diff.d" % str(case))
	t1 = data1[:,0]
	hz1 = data1[:,1]
	plot(hz_r_in,hz1)
	title('lissajous_absolute%s' % str(case))
	savefig("./out_lissajous/%slissajous12" % str(10000+case))
	# show()
