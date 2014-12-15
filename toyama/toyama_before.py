# coding:utf-8
#2011
import numpy as np
from pylab import *

data = np.genfromtxt("./NONAME03-01_edit.csv",delimiter=",",skiprows=1)

#ch1 X軸, ch1 Y軸, ch2 X軸, ch2 Y軸, ch5 移動距離, ch6 マーカー
l1 = data[:,0] #time[s]
l2 = data[:,1] #??
l3 = data[:,2] #自己比較X
l4 = data[:,3] #自己比較Y
l5 = data[:,4] #標準比較X
l6 = data[:,5] #標準比較Y
l7 = data[:,6] #ch5 max
l8 = data[:,7] #ch6 max

plot(l7,l8,label='marker')
plot(l7,l3,label='differential X')
plot(l7,l4,label='differential Y')
xlabel('distance')
ylabel('differential')
title('differential_toyama_before')
legend()
savefig('differentialXY')
show()

plot(l7,l8,label='marker')
plot(l7,l5,label='absolute X')
plot(l7,l6,label='absolute Y')
xlabel('distance')
ylabel('absolute')
title('absolute_toyama_before')
legend()
savefig('absoluteXY')
show()

