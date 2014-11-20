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

# aa = np.arange(1,101,1)
# bb = np.arange(1,101,1)
# AA,BB=np.meshgrid(aa,bb)

plt.xlabel('x')
plt.ylabel('z')


# im=plt.imshow(sig2,origin='lower',cmap=cm.RdYlBu) #blue red
# im=plt.imshow(sig2,origin='lower',cmap=cm.seismic) #seismic
# im=plt.imshow(sig2,origin='lower',cmap=cm.jet) #original
im=plt.imshow(sig2,origin='lower',cmap=cm.bwr) #blue red
# im=plt.imshow(sig2,origin='lower',cmap=cm.rainbow) #red_blue
# im=plt.imshow(sig2,origin='lower',cmap=cm.RdBu)

# plt.colorbar(im)
plt.savefig('sigmodel_xz.png')
# show()


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
data1 = np.loadtxt('./model_sig_xy.dat')
x1 = data1[:,0]
y1 = data1[:,1]
z1 = data1[:,2]
sig = data1[:,3]
sig2 = sig.reshape((101,101))

plt.xlabel('x')
plt.ylabel('y')

im=plt.imshow(sig2)

# plt.colorbar(im)
plt.savefig('sigmodel_xy.png')
# show()



# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
data1 = np.loadtxt('./model_myu_xz.dat')
x1 = data1[:,0]
y1 = data1[:,1]
z1 = data1[:,2]
sig = data1[:,3]
sig2 = sig.reshape((101,101))

plt.xlabel('x')
plt.ylabel('z')
plt.colorbar

im=plt.imshow(sig2,origin='lower')
1
plt.savefig('myumodel_xz.png')
# show()


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
data1 = np.loadtxt('./model_myu_xy.dat')
x1 = data1[:,0]
y1 = data1[:,1]
z1 = data1[:,2]
sig = data1[:,3]
sig2 = sig.reshape((101,101))


plt.xlabel('x')
plt.ylabel('y')

im=plt.imshow(sig2)

# plt.colorbar(im)
plt.savefig('myumodel_xy.png')
# show()