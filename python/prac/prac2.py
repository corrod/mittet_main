# coding:utf-8
#グローバル変数
x = 5


def addx(y):
    return x + y

print addx(10)


def setx(y):
    global x
    print x
    x = y
    print('x is %d' % x)
setx(10)

print x


def variable_args(*args, **kwargs):
    print 'args is ', args
    print 'kwargs is', kwargs

variable_args('ones', 'papapa', x=1, y=2, z=3)

va = variable_args
va('three', x=1, y=2)

# def funcname(params):
#     """concise one-line sentence describing the function.

#     Extended summary which can contain multiple paragraphs.
#     """
#     # function body
#     pass

# funcname?
