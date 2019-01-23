import math
from matplotlib import pyplot
import numpy as np


def plot(f, xmin, xmax, dx=0.001):
    x0 = xmin
    xplot = []
    fplot = []

    while x0 <= xmax:
        xplot.append(x0)
        fplot.append(f(x0))
        x0 += dx

    pyplot.plot(xplot, fplot)
    pyplot.show()

def cplot(x, f):
    fr = np.real(f)
    fi = np.imag(f)
    pyplot.plot(x, fr)
    pyplot.plot(x, fi)
    pyplot.show()
