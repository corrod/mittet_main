#coding:utf-8
#スライスとビューとは？
import numpy as np


#コピーとビュー
a = np.arange(10)
print 'a', a
#ビュー
b = a[::2]
b[0] = 12
print 'b', b
print 'a', a   # (!)

a = np.arange(10)
#コピー
b = a[::2].copy()  # force a copy
b[0] = 12
print 'b', b
print 'a', a

d = np.ones((100, 100))
d += d.T
print d

#次元の追加
z = np.array([1, 2, 3])
print z
print z[:, np.newaxis]
print z[np.newaxis, :]
