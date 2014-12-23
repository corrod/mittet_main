# coding:utf-8

from pylab import *
import numpy as np
import matplotlib.pyplot as plt


cases =range(19,103)

for i,case in enumerate(cases):
f =open("./out_residual/%sresidual_wave.d" % str(case),"w")
	data1 = np.loadtxt("../ver0030nocrack/%spattern2_23.d" % str(case))
	t1 = data1[:,0]
	hz1 = data1[:,1]
	data2 = np.loadtxt("./out_diff/%spattern2_23.d" % str(case))
	t2 = data2[:,0]
	hz2 = data2[:,1]
	residual = hz2 - hz1
	print >> f, t1, residual