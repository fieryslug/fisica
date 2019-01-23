import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

fig = plt.figure()
ax = plt.axes(xlim=(-8, 8), ylim=(-2, 2))
line, = ax.plot([], [], lw=2)
line2, = ax.plot([], [], lw=2)

tot_time = 20
intv = 20

def init():
    line.set_data([], [])
    line2.set_data([], [])
    return line,

def animate(i):
    x = np.linspace(-8, 8, 1000)
    y = np.sin(0.1*x*i)
    y2 = np.cos(0.1*x*i)
    line.set_data(x, y)
    line2.set_data(x, y2)
    return line, 

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=(tot_time * 1000) // intv, interval=intv, blit=True)
#anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()
