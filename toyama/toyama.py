# coding:utf-8
#2012
import numpy as np
from pylab import *
import matplotlib.pyplot as plt

data = np.genfromtxt("./ET_edit.csv",delimiter=",",skiprows=1)

#ch1 電流, ch2 電圧, ch 3GV, ch4 制圧機, ch5 回線速度, ch6 マーカー？?
l1 = data[:,0] #時間[s]
l2 = data[:,1] #ch1 max
l3 = data[:,2] #ch1 min
l4 = data[:,3] #ch2 max
l5 = data[:,4] #ch2 min
l6 = data[:,5] #ch3 max
l7 = data[:,6] #ch3 min
l8 = data[:,7] #ch4 max
l9 = data[:,8] #ch4 min
l10 = data[:,9] #ch5 max
l11 = data[:,10] #ch5 min
l12 = data[:,11] #ch6 max
l13 = data[:,12] #ch6 min

distance = l1 * l12
plot(l1,label="l1")
plot(l2,label="l2")
plot(l3,label="l3")
plot(l4,label="l4")
plot(l5,label="l5")
plot(l6,label="l6")
plot(l7,label="l7")
plot(l8,label="l8")
plot(l9,label="l9")
plot(l10,label="l10")
plot(l11,label="l11")
plot(l12,label="l12")
plot(l13,label="l13")

# plot(distance)
# plot(distance,l2)

legend()
savefig('toyama')
show()


# plot(l7,l3,label='differential X')
# plot(l7,l4,label='differential Y')
# xlabel('distance')
# ylabel('differential')

# legend()
# show()
# savefig('differential XY')


