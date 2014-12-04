
# coding:utf-8
from pylab import *
import numpy as np

import matplotlib.pyplot as plt


data1 = np.loadtxt('./model_sig_xz.dat')
x1 = data1[:,0]
y1 = data1[:,1]
z1 = data1[:,2]
sig = data1[:,3]
# sig2 = sig.reshape((101,101)) #101,101,101
sig2 = sig.reshape((71,101)) #101,61,71

# aa = np.arange(1,101,1)
# bb = np.arange(1,101,1)
# AA,BB=np.meshgrid(aa,bb)

fig=figure()
ax1=fig.add_subplot(111)

ax1.text('sea water')
ax1.set_xlabel('x')
ax1.set_ylabel('z')
# ax1.set_xticklabels([])
# ax1.set_yticklabels([])

# im=plt.imshow(sig2,origin='lower',cmap=cm.RdYlBu) #blue red
# im=plt.imshow(sig2,origin='lower',cmap=cm.seismic) #seismic
# im=plt.imshow(sig2,origin='lower',cmap=cm.jet) #original
plt.imshow(sig2,origin='lower',cmap=cm.bwr) #blue red
# im=plt.imshow(sig2,origin='lower',cmap=cm.rainbow) #red_blue
# im=plt.imshow(sig2,origin='lower',cmap=cm.RdBu)


plt.savefig('sigmodel_xz.png')
# show()







# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
data1 = np.loadtxt('./model_sig_xy.dat')
x1 = data1[:,0]
y1 = data1[:,1]
z1 = data1[:,2]
sig = data1[:,3]
# sig2 = sig.reshape((101,101))
sig2 = sig.reshape((61,101))


fig=figure()
ax1=fig.add_subplot(111)

ax1.set_title('model_sig_xy')

ax1.set_xlabel('x')
ax1.set_ylabel('y')
# ax1.set_xticklabels([])
# ax1.set_yticklabels([])

plt.imshow(sig2,origin='lower',cmap=cm.bwr) #blue red


plt.savefig('sigmodel_xy.png')
# show()



# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
data1 = np.loadtxt('./model_myu_xz.dat')
x1 = data1[:,0]
y1 = data1[:,1]
z1 = data1[:,2]
sig = data1[:,3]
# sig2 = sig.reshape((101,101))
sig2 = sig.reshape((71,101))

# plt.colorbar
fig=figure()
ax1=fig.add_subplot(111)

ax1.set_title('model_myu_xz')

ax1.set_xlabel('x')
ax1.set_ylabel('z')
# ax1.set_xticklabels([])
# ax1.set_yticklabels([])

plt.imshow(sig2,origin='lower',cmap=cm.bwr) #blue red

plt.savefig('myumodel_xz.png')
# show()


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
data1 = np.loadtxt('./model_myu_xy.dat')
x1 = data1[:,0]
y1 = data1[:,1]
z1 = data1[:,2]
sig = data1[:,3]
# sig2 = sig.reshape((101,101))
sig2 = sig.reshape((61,101))


fig=figure()
ax1=fig.add_subplot(111)

ax1.set_title('model_myu_xy')

ax1.set_xlabel('x')
ax1.set_ylabel('y')
# ax1.set_xticklabels([])
# ax1.set_yticklabels([])

plt.imshow(sig2,origin='lower',cmap=cm.bwr) #blue red


plt.savefig('myumodel_xy.png')
# show()