
# coding:utf-8
from pylab import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

cases = range(19,83)

for i, case in enumerate(cases):
	data1 = np.loadtxt("./out_diff/%spattern2_23diff.d" % str(case))
	t1 = data1[:,0]
	hz1 = data1[:,1]
	maxnumber=max(xrange(len(hz1)),key=lambda i:hz1[i])
	minnumber=min(xrange(len(hz1)),key=lambda i:hz1[i])
	# print case,maxnumber+1,max(hz1)
	print case,minnumber+1,min(hz1)



# data1 = np.loadtxt("./out_diff/43pattern2_23diff.d")
# t1 = data1[:,0]
# hz1 = data1[:,1]
# maxnumber=max(xrange(len(hz1)),key=lambda i:hz1[i])
# minnumber=min(xrange(len(hz1)),key=lambda i:hz1[i])
# print maxnumber+1,max(hz1), minnumber+1,min(hz1)


# data1 = np.loadtxt("./out_diff/45pattern2_23diff.d")
# t1 = data1[:,0]
# hz1 = data1[:,1]
# maxnumber=max(xrange(len(hz1)),key=lambda i:hz1[i])
# minnumber=min(xrange(len(hz1)),key=lambda i:hz1[i])
# print maxnumber+1,max(hz1), minnumber+1,min(hz1)


# data1 = np.loadtxt("./out_diff/47pattern2_23diff.d")
# t1 = data1[:,0]
# hz1 = data1[:,1]
# maxnumber=max(xrange(len(hz1)),key=lambda i:hz1[i])
# minnumber=min(xrange(len(hz1)),key=lambda i:hz1[i])
# print maxnumber+1,max(hz1), minnumber+1,min(hz1)


# data1 = np.loadtxt("./out_diff/49pattern2_23diff.d")
# t1 = data1[:,0]
# hz1 = data1[:,1]
# maxnumber=max(xrange(len(hz1)),key=lambda i:hz1[i])
# minnumber=min(xrange(len(hz1)),key=lambda i:hz1[i])
# print maxnumber+1,max(hz1), minnumber+1,min(hz1)


# data1 = np.loadtxt("./out_diff/51pattern2_23diff.d")
# t1 = data1[:,0]
# hz1 = data1[:,1]
# maxnumber=max(xrange(len(hz1)),key=lambda i:hz1[i])
# minnumber=min(xrange(len(hz1)),key=lambda i:hz1[i])
# print maxnumber+1,max(hz1), minnumber+1,min(hz1)
	

# data1 = np.loadtxt("./out_diff/53pattern2_23diff.d")
# t1 = data1[:,0]
# hz1 = data1[:,1]
# maxnumber=max(xrange(len(hz1)),key=lambda i:hz1[i])
# minnumber=min(xrange(len(hz1)),key=lambda i:hz1[i])
# print maxnumber+1,max(hz1), minnumber+1,min(hz1)
	

# data1 = np.loadtxt("./out_diff/55pattern2_23diff.d")
# t1 = data1[:,0]
# hz1 = data1[:,1]
# maxnumber=max(xrange(len(hz1)),key=lambda i:hz1[i])
# minnumber=min(xrange(len(hz1)),key=lambda i:hz1[i])
# print maxnumber+1,max(hz1), minnumber+1,min(hz1)
	

# data1 = np.loadtxt("./out_diff/57pattern2_23diff.d")
# t1 = data1[:,0]
# hz1 = data1[:,1]
# maxnumber=max(xrange(len(hz1)),key=lambda i:hz1[i])
# minnumber=min(xrange(len(hz1)),key=lambda i:hz1[i])
# print maxnumber+1,max(hz1), minnumber+1,min(hz1)
	

# data1 = np.loadtxt("./out_diff/59pattern2_23diff.d")
# t1 = data1[:,0]
# hz1 = data1[:,1]
# maxnumber=max(xrange(len(hz1)),key=lambda i:hz1[i])
# minnumber=min(xrange(len(hz1)),key=lambda i:hz1[i])
# print maxnumber+1,max(hz1), minnumber+1,min(hz1)
	

# data1 = np.loadtxt("./out_diff/61pattern2_23diff.d")
# t1 = data1[:,0]
# hz1 = data1[:,1]
# maxnumber=max(xrange(len(hz1)),key=lambda i:hz1[i])
# minnumber=min(xrange(len(hz1)),key=lambda i:hz1[i])
# print maxnumber+1,max(hz1), minnumber+1,min(hz1)
	

# data1 = np.loadtxt("./out_diff/63pattern2_23diff.d")
# t1 = data1[:,0]
# hz1 = data1[:,1]
# maxnumber=max(xrange(len(hz1)),key=lambda i:hz1[i])
# minnumber=min(xrange(len(hz1)),key=lambda i:hz1[i])
# print maxnumber+1,max(hz1), minnumber+1,min(hz1)
	

# data1 = np.loadtxt("./out_diff/65pattern2_23diff.d")
# t1 = data1[:,0]
# hz1 = data1[:,1]
# maxnumber=max(xrange(len(hz1)),key=lambda i:hz1[i])
# minnumber=min(xrange(len(hz1)),key=lambda i:hz1[i])
# print maxnumber+1,max(hz1), minnumber+1,min(hz1)
	

# data1 = np.loadtxt("./out_diff/67pattern2_23diff.d")
# t1 = data1[:,0]
# hz1 = data1[:,1]
# maxnumber=max(xrange(len(hz1)),key=lambda i:hz1[i])
# minnumber=min(xrange(len(hz1)),key=lambda i:hz1[i])
# print maxnumber+1,max(hz1), minnumber+1,min(hz1)
	

# data1 = np.loadtxt("./out_diff/69pattern2_23diff.d")
# t1 = data1[:,0]
# hz1 = data1[:,1]
# maxnumber=max(xrange(len(hz1)),key=lambda i:hz1[i])
# minnumber=min(xrange(len(hz1)),key=lambda i:hz1[i])
# print maxnumber+1,max(hz1), minnumber+1,min(hz1)
	

# data1 = np.loadtxt("./out_diff/71pattern2_23diff.d")
# t1 = data1[:,0]
# hz1 = data1[:,1]
# maxnumber=max(xrange(len(hz1)),key=lambda i:hz1[i])
# minnumber=min(xrange(len(hz1)),key=lambda i:hz1[i])
# print maxnumber+1,max(hz1), minnumber+1,min(hz1)
	

