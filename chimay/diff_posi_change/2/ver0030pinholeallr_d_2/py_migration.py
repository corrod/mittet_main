# coding:utf-8
import numpy as np
import scipy.signal
from pylab import *
from math import *
from cmath import *

nstep = 4351
# nstep = 4351

data_stac = np.loadtxt('./out_bp/19_2nd.dat')
cases = range(19 + 1, 113 + 1)
for i, case in enumerate(cases):
	data = np.loadtxt('./out_bp/%s_2nd.dat' % str(case))
	data_stac = data_stac + data

np.savetxt('migration.dat',data_stac)

fig = figure(figsize=(10,5))

c = data_stac
# c = np.transpose(data_stac)
a = len(c)
b = len(np.transpose(c))
xlabel('sp')
ylabel('time')
x = arange(1+18, b + 1+18)
y = arange(1, a + 1)
X,Y = meshgrid(x, y)
contourf(X, Y, c)
colorbar()
savefig('%smigration_space' % str(nstep))
show()