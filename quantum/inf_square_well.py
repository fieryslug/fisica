#Griffiths_P2_8

import math
from matplotlib import pyplot

a = 1
xmin = -a/4
xmax = 5*a/4

def eigenfunc(n, x):
    if 0 <= x <= a:
        return math.sqrt(2/a)*math.sin(n*math.pi*x/a)
    else:
        return 0

def eigenfuncfunc(n):
    def f(x):
        return eigenfunc(n, x)
    return f
        

def psi_0(x):
    if 0 <= x <= a/2:
        return math.sqrt(2/a)
    else:
        return 0

def add_func(f1, f2):
    def f(x):
        return f1(x) + f2(x)
    return f

def scal_mul_func(a, f1):
    def f(x):
        return a*f1(x)
    return f

def zero_func():
    def f(x):
        return 0
    return f

def in_prod(f1, f2):
    prod = 0
    x0 = xmin
    inc = 0.0001
    while x0 <= xmax:
        prod += f1(x0) * f2(x0) * inc
        x0 += inc
    return prod

def fourier(m, psi):

    f0 = zero_func()

    for j in range(1, m+1):
        eigenfuncj = eigenfuncfunc(j)
        cj = in_prod(eigenfuncj, psi)
        f0 = add_func(f0, scal_mul_func(cj, eigenfuncj))

    return f0
        
def plot(f, xmn=xmin, xmx=xmax):
    xplot = []
    psiplot = []

    x0 = xmin
    inc = 0.001
    while x0 <= xmax:
        xplot.append(x0)
        psiplot.append(f(x0))
        x0 += inc

    pyplot.plot(xplot, psiplot)
    pyplot.show()

