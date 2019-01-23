import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

fig = plt.figure()
ax = plt.axes(xlim=(-5, 5), ylim=(-3, 3))
liner, = ax.plot([], [], lw=0.5)
linei, = ax.plot([], [], lw=0.5)
linep, = ax.plot([], [], lw=0.7)

tot_time = 100
intv = 20

hb = 0.2
delta = 0.01
m = 100
p0 = 10

def init():
    liner.set_data([], [])
    linei.set_data([], [])
    linep.set_data([], [])
    return liner, linei, linep

def animate(i):

    x = np.linspace(-5, 5, 5000)
    t = i*intv/1000
    
    y1 = (np.pi**(1/2)*(delta+1j*hb*t/(m*delta))) ** -1/2
    y2 = np.exp(-(x-p0*t/m)**2/(2*delta**2*(1+1j*hb*t/(m*delta**2))))
    y3 = np.exp((1j*p0/hb)*(x-p0*t/(2*m)))
    y = y1 * y2 * y3 * 2
    yr = np.real(y)
    yi = np.imag(y)
    liner.set_data(x, yr)
    linei.set_data(x, yi)
    linep.set_data(x, y * np.conj(y))
    return linei, liner, linep

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=(tot_time * 1000) // intv, interval=intv, blit=True)
#anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()
