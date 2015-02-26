# coding:utf-8
import numpy as np
import scipy.signal
from pylab import *
from math import *
from cmath import *
import scipy.interpolate

data = np.loadtxt('test.dat')
#x_t =linspace(5,700,15)
x = data[:,0]
real = data[:,1]
imag = data[:,2]
amp = (real[:]**2 + imag[:]**2)**0.5
phase = arctan(imag/real) /3.14*180


# sp1 = scipy.interpolate.InterpolatedUnivariateSpline(x,amp)
# sx1 = np.linspace(5,150,50)
# sy1 = sp1(sx1)

plot(x[0:11],amp[0:11],'o',label='theoretical')
plot(x[0:12],amp[0:12]+0.00000001,label='simulation',lw=1)
# plot(sx1,sy1)
xlabel('Distance [m]',fontsize=20)
ylabel('Amplitude [A/m]',fontsize=20)
yscale('log')
title('Amplitude')
legend()
savefig('amplitude_test')
show()


sp = scipy.interpolate.InterpolatedUnivariateSpline(x,phase)
sx = np.linspace(5,150,100)
sy=sp(sx)
plot(x[0:11],phase[0:11],'o',label='theoretical')
# plot(x[0:12],phase[0:12]+0.2,label='simulation',lw=2)
plot(sx,sy+0.3)
xlabel('Distance [m]',fontsize=20)
ylabel('Phase',fontsize=20)
title('Phase')
legend()
savefig('phase_test')
show()
