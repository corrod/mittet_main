# coding:utf-8
import numpy as np
import scipy.signal
from pylab import *
from math import *
from cmath import *

# x = arange(0,2*pi,2*pi/360)
x = linspace(0,2*pi,360)

a = np.cos(x)
d = real(scipy.signal.hilbert(a,axis=-1))
e = imag(scipy.signal.hilbert(a,axis=-1))

plot(x,a,label='cos')
plot(x,d,label='hilbert_real')
plot(x,e,label='hilbert_imag')
title('hilbert transform cos')
grid(True)
legend()
# savefig('hilbert_cos.png')
show()


a = np.sin(x)
d = real(scipy.signal.hilbert(a,axis=-1))
e = imag(scipy.signal.hilbert(a,axis=-1))

plot(x,a,label='sin')
plot(x,d,label='hilbert_real')
plot(x,e,label='hilbert_imag')
title('hilbert transform sin')
grid(True)
legend()
# savefig('hilbert_sin.png')
show()