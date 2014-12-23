import numpy as np

g = range(1,86)
np.savetxt('pycordinate.dat',g)

cases = range(19,103)

data2 = np.loadtxt("./out_residual/19residual_wave_abso.d")
hz2 = data2[:,1]
for i, case in enumerate(cases):
	data1 = np.loadtxt("./out_residual/%sresidual_wave_abso.d" % str(case))
	t1 = data1[:,0]
	hz1 = data1[:,1]
	hz2 = np.c_[hz2,hz1]
np.savetxt('combine_residual_abso.dat',hz2)
# f = open('combine_residual_abso.dat')
# print >> f,hz2
# f.close()


data22 = np.loadtxt("./out_residual/19residual_wave_diff.d")
hz22 = data22[:,1]
for i, case in enumerate(cases):
	data11 = np.loadtxt("./out_residual/%sresidual_wave_diff.d" % str(case))
	t11 = data11[:,0]
	hz11 = data11[:,1]
	hz22 = np.c_[hz22,hz11]
np.savetxt('combine_residual_diff.dat',hz22)
# f = open('combine_residual_abso.dat')
# print >> f,hz2
# f.close()
