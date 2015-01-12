# coding:utf-8

import numpy as np
from pylab import *

inputdata = np.loadtxt('./signal.d')
t_in = inputdata[:,0]
hz_r_in = inputdata[:,1]
hz_i_in = inputdata[:,2]

cases = range(19,114)

for i, case in enumerate(cases):
	data1 = np.loadtxt("./out_diff/%spattern2_12diff.d" % str(case))
	t1 = data1[:,0]
	hz1 = data1[:,1]
	plot(hz_r_in,hz1)
	title('lissajous_absolute%s' % str(case))
	savefig("./out_lissajous/%slissajous12" % str(10000+case))
	# show()

for i, case in enumerate(cases):
	data2 = np.loadtxt("./out_diff/%spattern2_23diff.d" % str(case))
	t2 = data2[:,0]
	hz2 = data2[:,1]
	plot(hz_r_in,hz2)
	title('lissajous_differential%s' % str(case) )
	# xlim(xmax=7000000)
	savefig("./out_lissajous/%slissajous23" % str(10000+case))
	# show()