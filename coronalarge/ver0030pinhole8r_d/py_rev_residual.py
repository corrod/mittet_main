# coding:utf-8
# residual_wave_abso.d residal_wave_diff.d wo reverse time
import numpy as np
import matplotlib.pyplot as plt

cases = range(19,114)

#reverse
for i, case in enumerate(cases):
	data = np.loadtxt("./out_residual/%sresidual_wave_abso.d" % str(case))
	t = data[:,0]
	resi_data = data[:,1]
	rev_resi_data = resi_data[::-1] #reverse [::-1]
	rev_resi_data = np.c_[t,rev_resi_data]
	np.savetxt("./out_residual_rev/%sresidual_wave_abso_rev.d" % str(case),rev_resi_data)

for i, case in enumerate(cases):
	data = np.loadtxt("./out_residual/%sresidual_wave_diff.d" % str(case))
	t = data[:,0]
	resi_data = data[:,1]
	rev_resi_data = resi_data[::-1]
	rev_resi_data = np.c_[t,rev_resi_data]
	np.savetxt("./out_residual_rev/%sresidual_wave_diff_rev.d" % str(case),rev_resi_data)

#combine
data22 = np.loadtxt("./out_residual_rev/19residual_wave_abso_rev.d")
hz22 = data22[:,1]
for i, case in enumerate(cases):
	data11 = np.loadtxt("./out_residual_rev/%sresidual_wave_abso_rev.d" % str(case))
	t = data11[:,0]
	hz11 = data11[:,1]
	hz22 = np.c_[hz22,hz11]
	np.savetxt("combine_residual_abso_rev.dat",hz22)

data22 = np.loadtxt("./out_residual_rev/19residual_wave_diff_rev.d")
hz22 = data22[:,1]
for i, case in enumerate(cases):
	data11 = np.loadtxt("./out_residual_rev/%sresidual_wave_diff_rev.d" % str(case))
	t = data11[:,0]
	hz11 = data11[:,1]
	hz22 = np.c_[hz22,hz11]
	np.savetxt("combine_residual_diff_rev.dat",hz22)
