from pylab import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

data1 = np.loadtxt('./model_sig_xz.dat')
x1 = data1[:,0]
y1 = data1[:,1]
z1 = data1[:,2]
sig = data1[:,3]
sig2 = sig.reshape((101,101))

aa = np.arange(1,101,1)
bb = np.arange(1,101,1)
AA,BB=np.meshgrid(aa,bb)

plt.xlabel('x')
plt.ylabel('z')

plt.imshow(sig2)
# plt.contourf(AA,BB,sig2,2,alpha=.75,cmap=plt.cm.hot)
plt.colorbar
plt.savefig('sigmodel_xz.png')
# show()

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
data1 = np.loadtxt('./model_sig_xy.dat')
x1 = data1[:,0]
y1 = data1[:,1]
z1 = data1[:,2]
sig = data1[:,3]
sig2 = sig.reshape((101,101))

aa = np.arange(1,101,1)
bb = np.arange(1,101,1)
AA,BB=np.meshgrid(aa,bb)

plt.xlabel('x')
plt.ylabel('z')

plt.imshow(sig2)

plt.colorbar
plt.savefig('sigmodel_xy.png')
# show()

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
data1 = np.loadtxt('./model_myu_xy.dat')
x1 = data1[:,0]
y1 = data1[:,1]
z1 = data1[:,2]
sig = data1[:,3]
sig2 = sig.reshape((101,101))

aa = np.arange(1,101,1)
bb = np.arange(1,101,1)
AA,BB=np.meshgrid(aa,bb)

plt.xlabel('x')
plt.ylabel('z')

plt.imshow(sig2)

plt.colorbar
plt.savefig('myumodel_xy.png')
# show()


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
data1 = np.loadtxt('./model_myu_xz.dat')
x1 = data1[:,0]
y1 = data1[:,1]
z1 = data1[:,2]
sig = data1[:,3]
sig2 = sig.reshape((101,101))

aa = np.arange(1,101,1)
bb = np.arange(1,101,1)
AA,BB=np.meshgrid(aa,bb)

plt.xlabel('x')
plt.ylabel('z')
plt.colorbar

plt.imshow(sig2)
# plt.hot()
# plt.gray()
plt.savefig('myumodel_xz.png')
# show()