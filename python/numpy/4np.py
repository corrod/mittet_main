#coding utf-8

import numpy as np
from matplotlib import pyplot as plt

data = np.loadtxt('populations.txt')
print data
year, hares, lynxes, carrots = data.T
plt.axes([0.2, 0.1, 0.5, 0.8])
plt.plot(year, hares, year, lynxes, year, carrots)
plt.legend(('Hare', 'Lynx', 'Carrot'), loc=(1.05, 0.5))
plt.show()

populations = data[:, 1:]
print populations
print 'mean', populations.mean(axis=0)
print 'sed', populations.std(axis=0)
print 'argmax', np.argmax(populations, axis=1)
