#coding utf-8

import numpy as np
from matplotlib import pyplot as plt

from scipy.integrate import quad

res, err = quad(np.sin, 0, np.pi/2)
np.allclose(res, 1)
np.allclose(err, 1 - res)


def calc_derivative(ypos, time, counter_arr):
    counter_arr += 1
    return -2 * ypos

counter = np.zeros((1,), dtype=np.uint16)

from scipy.integrate import odeint
time_vec = np.linspace(0, 4, 40)
yvec, info = odeint(calc_derivative, 1, time_vec,
                    args=(counter,), full_output=True)

print counter
print info['nfe'][:10]
