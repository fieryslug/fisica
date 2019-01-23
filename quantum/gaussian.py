import plotter as pltr
from matplotlib import pyplot as plt
import numpy as np

delta = 0.001
a = 0
p = 10
hb = 1
k0 = p / hb

def gauss(x):
    return (np.pi*delta**2)**(-1/4) * np.exp(1j*k0*(x+a)) * np.exp(-(x+a)**2/(2*delta**2))


x = np.linspace(-2, 2, 2000)
y = gauss(x)
yr = np.real(y)
yi = np.imag(y)

plt.plot(x, yr, x, yi)
plt.show()
