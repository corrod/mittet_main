# coding:utf-8

import numpy as np
import scipy as scipy
print "hello world"
a = 3
b = 2*a
print "type(b)", type(b)
print b
a * b
b = "hello world"
print type(b)
b + b
2 * b
print b + b
a = 1
if 2**2 == 4:
    print 'obvious!'
if a == 1:
    print(1)
elif a == 2:
    print(2)
else:
    print('a lot')
print 3

for i in range(4):
    print(i)

# for word in ('cool', 'powerful', 'readable'):
#     print('python is %s' % word)

# for word in ('cool', 'powerful', 'readable'):
#     print('python is %s' % word)


z = 1 + 1j
while abs(z) < 100:
    z = z**2 + 1
print(z)

z = 1 + 1j
while abs(z) < 100:
    if z.imag == 0:
        break
    z = z**2 + 1
    print z

a = [1, 0, 2, 4]
for element in a:
    if element == 0:
        continue
    print 1. / element

# 反復回数の追跡
words = ('cool', 'powerful', 'readable')
for index, item in enumerate(words):
    print index, item

# # wollis の公式
# for i in range(100):
#     pi = 2. * (4. * i**2.) / (4. * i**2. - 1.)
#     print pi


def test():
    print('in test function')

test()


def disk_area(radius):
    return 3.14 * radius * radius

print disk_area(1.5)


foo = [1, 2, 3]
print [2 * x for x in foo]


poo = [1, 2, 3, 4, 5, 6]
print [x for x in poo if x % 3 == 0]


for k in xrange(10):
    print k


for h in range(10):
    print h

g = xrange(19)
print g

for kk in g:
    print kk


# print dir(poo)


# x = ('foo', 'bar')
# print x[0]
# x[0] = 'hoge'

def hoo(x=[]):
    x.append(1)
    print x

hoo()
hoo()
hoo()


def hoo(x=[]):
    if not x:
        x = []
    x.append(1)
    print x

hoo()
hoo()
hoo()


class qoo():
    @staticmethod
    def x():
        return 3


def too(*args):
    print args

too(4, 5, 6)


def koo(**kwargs):
    print kwargs

koo(x=3, y=0)


def paa(x, y=0, z=0):
    print x, y, z
paa(1, 2, 4)
paa(1, 2)
paa(1)
paa(0)
paa(z=3, y=2, x=1)


def try_to_modify(v, w, z):
    v = 23
    w.append(42)
    z = [99]
    print(v)
    print(w)
    print(z)

a = 77
b = [99]
c = [28]

try_to_modify(a, b, c)

print a, b, c
