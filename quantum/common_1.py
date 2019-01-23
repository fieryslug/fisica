import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from time import time

hbar = 1

class Particle:
    def __init__(self, mass, xmin, xmax, pres=1000):
        self.x = np.linspace(xmin, xmax, pres)
        self.mass = mass
        self.xmin = xmin
        self.xmax = xmax
        self.dx = (xmax - xmin) / (pres-1)
        self.len = pres
        self.state = 0 * self.x
        self.time = 0
        self.poten = 0 * np.ndarray(self.len)
        
    def Hamilton(self):
        return -(hbar**2/(2*self.mass)) * self.ddx(self.ddx(self.state)) + self.poten * self.state

    def ddx(self, func):
        
        return np.around(np.gradient(func, edge_order=2) / self.dx, 10)

    def ddx_deprct(self, func):
        der = 1j*np.ndarray(self.len)
        der[0] = (func[1] - func[0]) / self.dx
        der[-1] = (func[-1] - func[-2]) / self.dx
        for i in range(1, self.len-1):
            dfdx = (func[i+1] - func[i-1]) / (2*self.dx)
            if np.isnan(dfdx):
                print('i:: {} - {} = {}'.format(func[i+1], func[i-1], func[i+1]-func[i-1]))

        return np.around(der, 5)
        
    def update(self, dt):
        self.state += -(1j/hbar) * self.Hamilton() * dt
        self.time += dt

def gauss_packet(x, p0=0, delta=1):
    return np.exp(1j*p0*x/hbar) * np.exp(-x**2/(2*delta**2)) / ((np.pi*delta**2)**(1/4))

def animator(ptc, dt, tot_time):
    fig = plt.figure()
    ax = plt.axes(xlim=(-8, 8), ylim=(-3, 3))
    liner, = ax.plot([], [], lw=0.5)
    linei, = ax.plot([], [], lw=0.5)
    linep, = ax.plot([], [], lw=0.7)

    

    def init():
        liner.set_data([], [])
        linei.set_data([], [])
        linep.set_data([], [])
        return liner, linei, linep,

    def animate(i):
        #print(i)
        ptc.update(dt)
        stater = np.real(ptc.state)
        statei = np.imag(ptc.state)
        p = stater**2 + statei**2
        liner.set_data(ptc.x, stater)
        linei.set_data(ptc.x, statei)
        linep.set_data(ptc.x, p)
        
        return linep,

    t0 = time()
    animate(0)
    t1 = time()
    intv = int(1000 * dt - (t1-t0))
    print(intv)

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=(tot_time*1000)//intv, interval=intv, blit=True)
    plt.show()

