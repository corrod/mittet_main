#coding utf-8

import numpy as np
from matplotlib import pyplot as plt
from scipy import signal

t = np.linspace(0, 5, 100)
x = t + np.random.normal(size=100)

plt.plot(t, x, linewidth=3)
plt.plot(t, signal.detrend(x), linewidth=3)
plt.show()

x = np.sin(t)
plt.plot(t, x, linewidth=3)
plt.plot(t[::2], signal.resample(x, 50), 'ko')
plt.show()
