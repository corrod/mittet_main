# coding:utf-8

from pylab import *
import numpy as np
import matplotlib.pyplot as plt

cwa = 1767.77
dt = 1.258468585615218E-006

h1 = dt*8.* cwa
h2 = dt*16.* cwa
print 'wall thinning h1, h2', h1, h2
print 'wall thinning true h1, h2 ', 4.5E-003, 4.5E-03*2.
c1 = 4.5E-003 / (dt*8.)
c2 = (4.5E-003*2.) / (dt*16.)
print 'velocity c1, c2', c1, c2