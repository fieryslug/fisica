import math
from matplotlib import pyplot


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
