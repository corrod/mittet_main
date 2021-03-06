#coding utf-8

import numpy as np
from matplotlib import pyplot as plt

# p = np.polynomial.Polynomial([-1, 2, 3])  # coefs in different order!
# print p

x = np.linspace(-1, 1, 2000)
y = np.cos(x) + 0.3*np.random.rand(2000)
p = np.polynomial.Chebyshev.fit(x, y, 90)
t = np.linspace(-1, 1, 200)
plt.plot(x, y, 'r.')
plt.plot(t, p(t), 'k-', lw=3)
plt.show()
