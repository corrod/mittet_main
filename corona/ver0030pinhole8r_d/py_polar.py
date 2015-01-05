# coding:utf-8
from pylab import *
import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt('./out_diff/45pattern2_23diff.d')
x1 = data1[:,0]
y1 = data1[:,1]

#theta = omega *t
theta = 2.* np.pi* 2000. *x1
theta2 = theta -10

#radian = do * pi/ 180
radian = theta * np.pi/180

ax = plt.axes(polar=True)
# ax.scatter(theta,y1)
ax.plot(theta,y1,'o',ls='-',ms=3,markevery=6)
ylim(ymax=1e10)
plt.xlabel('45diff23')
plt.savefig('polar45_23.png')
# plt.show()


# data2 = np.loadtxt('./out_diff/45pattern2_23diff.d')
# x2 = data2[:,0]
# y2 = data2[:,1]

# #theta = omega *t
# theta2 = 2.* np.pi* 2000. *x2 -10

# #radian = do * pi/ 180
# radian = theta * np.pi/180

# ax = plt.axes(polar=True)
# ax.scatter(theta,y2)

# plt.xlabel('45diff23kaiten')
# plt.savefig('polar45_23kaiten.png')
# # plt.show()


# data3 = np.loadtxt('./out_diff/45pattern2_23diff.d')
# x3 = data1[:,0]
# y3 = data1[:,1]

# #theta = omega *t
# theta3 = 2.* np.pi* 2000. *x3 -20

# #radian = do * pi/ 180
# radian = theta * np.pi/180

# ax = plt.axes(polar=True)
# ax.scatter(theta2,y1)

# plt.xlabel('45diff23kaiten2')
# plt.savefig('polar45_23kaiten2.png')
# # plt.show()




# data1 = np.loadtxt('./out_diff/45pattern2_12diff.d')
# x1 = data1[:,0]
# y1 = data1[:,1]
# #theta = omega *t
# theta = 2.* np.pi* 2000. *x1
# #radian = do * pi/ 180
# radian = theta * np.pi/180
# ax = plt.axes(polar=True)

# ax.scatter(theta,y1)
# plt.xlabel('45diff12')
# plt.savefig('polar45_12.png')
# # plt.show()