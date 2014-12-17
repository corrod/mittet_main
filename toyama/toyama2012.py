# coding:utf-8
#2012
import numpy as np
from pylab import *
import matplotlib.pyplot as plt

data = np.genfromtxt("./hakei2012/AUTO0001.CSV",delimiter=",",skiprows=1)

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


plot(l1,label="l1 time")
legend()
show()

plot(l1*l10,label='l1*l10')
legend()
show()

plot(l1*l10,l2,label='l1*l10,l2')
legend()
show()

plot(l1*l10,l4,label='l1*l10,l4')
legend()
show()


plot(l1,l12,label='marker')
plot(l1,l4,label='voltage')
legend()
savefig('marker-voltage')
show()

plot(l1,l12,label='marker')
plot(l1,l2,label='electric current')
legend()
savefig('marker-electric current')
show()

plot(l1,l12,label='marker')
plot(l1,l2,label='electric current')
plot(l1,l4,label='voltage')
legend()
savefig('electriccurrent-voltage')
show()



fig=figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ax1.plot(l2,label="l2")
ax2.plot(l3,label="l3")
ax1.legend()
ax2.legend()
title('electric current')
savefig('electric current')
show()

fig=figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ax1.plot(l4,label="l4")
ax2.plot(l5,label="l5")
ax1.legend()
ax2.legend()
title('voltage')
savefig('voltage')
show()

fig=figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ax1.plot(l6,label="l6")
ax2.plot(l7,label="l7")
ax1.legend()
ax2.legend()
title('GV')
savefig('GV')
show()


fig=figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ax1.plot(l8,label="l8")
ax2.plot(l9,label="l9")
ax1.legend()
ax2.legend()
title('seiatsuki??')
savefig('seiatsuki??')
show()

fig=figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ax1.plot(l10,label="l10")
ax2.plot(l11,label="l11")
ax1.legend()
ax2.legend()
title('rotational speed')
savefig('rotational speed')
show()

fig=figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ax1.plot(l12,label="l12") #marker max
ax2.plot(l13,label="l13") #marker min
ax1.legend()
ax2.legend()
title('marker')
savefig('marker')
show()



# plot(l7,l3,label='differential X')
# plot(l7,l4,label='differential Y')
# xlabel('distance')
# ylabel('differential')

# legend()
# show()
# savefig('differential XY')


