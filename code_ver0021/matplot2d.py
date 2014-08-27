import pylab as pl
import numpy as np

i, j, k, sig = pl.loadtxt('model_sig.dat', unpack=True)
# i, j, k, myu = pl.loadtxt('model_myu.dat')


# def f(x, y):
#     return (1 - x / 2 + x**5 + y**3) * np.exp(-x**2 -y**2)


# n = 256
# x = np.linspace(-3, 3, n)
# y = np.linspace(-3, 3, n)
# X, Y = np.meshgrid(x, y)

I, K = np.meshgrid(i, k)

pl.axes([0.025, 0.025, 0.95, 0.95])

pl.contourf(I, K, sig, 8, alpha=.75, cmap=pl.cm.hot)
C = pl.contour(I, K, sig, 8, colors='black', linewidth=.5)
# pl.contourf(X, Y, f(X, Y), 8, alpha=.75, cmap=pl.cm.hot)
# C = pl.contour(X, Y, f(X, Y), 8, colors='black', linewidth=.5)
pl.clabel(C, inline=1, fontsize=10)

pl.xticks(())
pl.yticks(())
pl.show()


# pylab.subplot(211)
