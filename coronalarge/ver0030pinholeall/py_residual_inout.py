# coding:utf-8

from pylab import *
import numpy as np
import matplotlib.pyplot as plt

cases =range(19,104)

# residual differential
for i,case in enumerate(cases):
	f =open("./out_residual/%sresidual_wave_diff.d" % str(case),"w")
	data1 = np.loadtxt("../ver0030nocrack/out_diff/%spattern2_23diff.d" % str(case))
	t1 = data1[:,0]
	hz1 = data1[:,1]

	data2 = np.loadtxt("./out_diff/%spattern2_23diff.d" % str(case))
	t2 = data2[:,0]
	hz2 = data2[:,1]
	residual = hz2 - hz1

	for ix in range(len(residual)):
		print >> f, t1[ix], residual[ix]
	f.close

# f = open("./coodinate.txt","w")
# for  case in enumerate(cases):
# 	print >> f, case
# f.close


# residual absolute
for i,case in enumerate(cases):
	f =open("./out_residual/%sresidual_wave_abso.d" % str(case),"w")
	data1 = np.loadtxt("../ver0030nocrack/out_diff/%spattern2_12diff.d" % str(case))
	t1 = data1[:,0]
	hz1 = data1[:,1]

	data2 = np.loadtxt("./out_diff/%spattern2_12diff.d" % str(case))
	t2 = data2[:,0]
	hz2 = data2[:,1]
	residual = hz2 - hz1

	for ix in range(len(residual)):
		print >> f, t1[ix], residual[ix]
	f.close

	# plot(t1,residual)
	# show()


#subplot residual differential
f = figure()
subplots_adjust(hspace=0.001)

ax1 = subplot(911)
ax1.set_title('residual_differential')
data1 = np.loadtxt('./out_residual/19residual_wave_diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
# ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('19')

ax1 = subplot(912)
data1 = np.loadtxt('./out_residual/29residual_wave_diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('29')

ax1 = subplot(913)
data1 = np.loadtxt('./out_residual/39residual_wave_diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('39')

ax1 = subplot(914)
data1 = np.loadtxt('./out_residual/49residual_wave_diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('49')

ax1 = subplot(915)
data1 = np.loadtxt('./out_residual/59residual_wave_diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('59')

ax1 = subplot(916)
data1 = np.loadtxt('./out_residual/69residual_wave_diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('69')

ax1 = subplot(917)
data1 = np.loadtxt('./out_residual/79residual_wave_diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('79')

ax1 = subplot(918)
data1 = np.loadtxt('./out_residual/89residual_wave_diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('89')

ax1 = subplot(919)
data1 = np.loadtxt('./out_residual/99residual_wave_diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
ax1.set_yticklabels([])
# ylim(ymax=2e11)
ax1.set_ylabel('99')

plt.xlabel('time [s]')

savefig('residual_differential')
# show()



#### differential
f = figure()
subplots_adjust(hspace=0.001)

ax1 = subplot(911)
ax1.set_title('differential')
data1 = np.loadtxt('./out_diff/19pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
# ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('19')

ax1 = subplot(912)
data1 = np.loadtxt('./out_diff/29pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('29')

ax1 = subplot(913)
data1 = np.loadtxt('./out_diff/39pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('39')

ax1 = subplot(914)
data1 = np.loadtxt('./out_diff/49pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('49')

ax1 = subplot(915)
data1 = np.loadtxt('./out_diff/59pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('59')

ax1 = subplot(916)
data1 = np.loadtxt('./out_diff/69pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('69')

ax1 = subplot(917)
data1 = np.loadtxt('./out_diff/79pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('79')

ax1 = subplot(918)
data1 = np.loadtxt('./out_diff/89pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('89')

ax1 = subplot(919)
data1 = np.loadtxt('./out_diff/99pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
ax1.set_yticklabels([])
# ylim(ymax=2e11)
ax1.set_ylabel('99')

plt.xlabel('time [s]')

savefig('differential')
# show()


#residual absolute
f = figure()
subplots_adjust(hspace=0.001)

ax1 = subplot(911)
ax1.set_title('residual_absolute')
data1 = np.loadtxt('./out_residual/19residual_wave_abso.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
# ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('19')

ax1 = subplot(912)
data1 = np.loadtxt('./out_residual/29residual_wave_abso.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('29')

ax1 = subplot(913)
data1 = np.loadtxt('./out_residual/39residual_wave_abso.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('39')

ax1 = subplot(914)
data1 = np.loadtxt('./out_residual/49residual_wave_abso.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('49')

ax1 = subplot(915)
data1 = np.loadtxt('./out_residual/59residual_wave_abso.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('59')

ax1 = subplot(916)
data1 = np.loadtxt('./out_residual/69residual_wave_abso.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('69')

ax1 = subplot(917)
data1 = np.loadtxt('./out_residual/79residual_wave_abso.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('79')

ax1 = subplot(918)
data1 = np.loadtxt('./out_residual/89residual_wave_abso.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('89')

ax1 = subplot(919)
data1 = np.loadtxt('./out_residual/99residual_wave_abso.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
ax1.set_yticklabels([])
# ylim(ymax=2e11)
ax1.set_ylabel('99')

plt.xlabel('time [s]')

savefig('residual_absolute')
# show()


#### absolute
f = figure()
subplots_adjust(hspace=0.001)

ax1 = subplot(911)
ax1.set_title('absolute')
data1 = np.loadtxt('./out_diff/19pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
# ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('19')

ax1 = subplot(912)
data1 = np.loadtxt('./out_diff/29pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('29')

ax1 = subplot(913)
data1 = np.loadtxt('./out_diff/39pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('39')

ax1 = subplot(914)
data1 = np.loadtxt('./out_diff/49pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('49')

ax1 = subplot(915)
data1 = np.loadtxt('./out_diff/59pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('59')

ax1 = subplot(916)
data1 = np.loadtxt('./out_diff/69pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('69')

ax1 = subplot(917)
data1 = np.loadtxt('./out_diff/79pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('79')

ax1 = subplot(918)
data1 = np.loadtxt('./out_diff/89pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('89')

ax1 = subplot(919)
data1 = np.loadtxt('./out_diff/99pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
ax1.set_yticklabels([])
# ylim(ymax=2e11)
ax1.set_ylabel('99')

plt.xlabel('time [s]')

savefig('absolute')
# show()

#### absolute nocrack
f = figure()
subplots_adjust(hspace=0.001)

ax1 = subplot(911)
ax1.set_title('absolute nocrack')
data1 = np.loadtxt('../ver0030nocrack/out_diff/19pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
# ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('19')

ax1 = subplot(912)
data1 = np.loadtxt('../ver0030nocrack/out_diff/29pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('29')

ax1 = subplot(913)
data1 = np.loadtxt('../ver0030nocrack/out_diff/39pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('39')

ax1 = subplot(914)
data1 = np.loadtxt('../ver0030nocrack/out_diff/49pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('49')

ax1 = subplot(915)
data1 = np.loadtxt('../ver0030nocrack/out_diff/59pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('59')

ax1 = subplot(916)
data1 = np.loadtxt('../ver0030nocrack/out_diff/69pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('69')

ax1 = subplot(917)
data1 = np.loadtxt('../ver0030nocrack/out_diff/79pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('79')

ax1 = subplot(918)
data1 = np.loadtxt('../ver0030nocrack/out_diff/89pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
# ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('89')

ax1 = subplot(919)
data1 = np.loadtxt('../ver0030nocrack/out_diff/99pattern2_12diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
ax1.set_yticklabels([])
# ylim(ymax=2e11)
ax1.set_ylabel('99')

plt.xlabel('time [s]')

savefig('absolute_nocrack')
# show()


#### differential nocrack
f = figure()
subplots_adjust(hspace=0.001)

ax1 = subplot(911)
ax1.set_title('differential_nocrack')
data1 = np.loadtxt('../ver0030nocrack/out_diff/19pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
ylim(ymax=2e11)
# ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('19')

ax1 = subplot(912)
data1 = np.loadtxt('../ver0030nocrack/out_diff/29pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('29')

ax1 = subplot(913)
data1 = np.loadtxt('../ver0030nocrack/out_diff/39pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('39')

ax1 = subplot(914)
data1 = np.loadtxt('../ver0030nocrack/out_diff/49pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('49')

ax1 = subplot(915)
data1 = np.loadtxt('../ver0030nocrack/out_diff/59pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('59')

ax1 = subplot(916)
data1 = np.loadtxt('../ver0030nocrack/out_diff/69pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('69')

ax1 = subplot(917)
data1 = np.loadtxt('../ver0030nocrack/out_diff/79pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('79')

ax1 = subplot(918)
data1 = np.loadtxt('../ver0030nocrack/out_diff/89pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
ylim(ymax=2e11)
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel('89')

ax1 = subplot(919)
data1 = np.loadtxt('../ver0030nocrack/out_diff/99pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]
plt.plot(x1,y1)
ax1.set_yticklabels([])
ylim(ymax=2e11)
ax1.set_ylabel('99')

plt.xlabel('time [s]')

savefig('differential_nocrack')
# show()