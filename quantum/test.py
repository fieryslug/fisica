import numpy as np
import plotter as pltr
from scipy.fftpack import fft, ifft
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import common as mq
import plotter as pltr
import time


fig = plt.figure()
ax = fig.add_subplot()


ax = plt.axes(xlim = (-15, 15), ylim = (-2, 2))
liner, = ax.plot([], [], lw=1)
linei, = ax.plot([], [], lw=1)
linep, = ax.plot([], [], lw=1)
linev, = ax.plot([], [], lw=1)

x = np.linspace(-10, 10, 2**12)
psi_x0 = mq.gauss_packet(x, p0=15)
psi_x0[:100] = 0
psi_x0[-100:] = 0
ptc = mq.ParticleInABox(x, psi_x0, 5, dt=0.01)

def init():
    liner.set_data([], [])
    linei.set_data([], [])
    linep.set_data([], [])
    return linep,

def animate(i):
    ptc.update()
    linep.set_data(ptc.xax, ptc.prob)
    linev.set_data(ptc.xax, ptc.poten/50)
    print(ptc.time)
    return linep, linev,

t0 = time.time()
animate(0)
t1 = time.time()

intv = 1000 * ptc.dt - (t1-t0)

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=300, interval=intv, blit=True)
plt.show()
