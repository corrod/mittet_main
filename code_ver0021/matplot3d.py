from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

x = np.arange(1, 50, 0.25)
y = np.arange(1, 50, 0.25)
z = np.sin(X) + np.cos(Y)
X, Y = np.meshgrid(x, y)
fig = plt.figure()
ax = Axes3D(fig)
ax.plot_wireframe(X, Y, Z)

plt.show()
