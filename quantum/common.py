import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import plotter as pltr
from scipy.fftpack import fft, ifft, fftshift
import time

hbar = 1.0

def gauss_packet(x, p0=0, delta=1):
    return np.exp(1j*p0*x/hbar) * np.exp(-x**2/(2*delta**2)) / ((np.pi*delta**2)**(1/4))

def plane_wave(x, p0):
    return np.exp(1j*p0*x/hbar)

def wall(width, x):
    y = 0 * x
    N = len(y)
    for i in range(width):
        y[i] = 1e8
        y[N-1-i] = 1e8
    return y

def animator(ptc, tot_time, drawPotent=True): 
    fig = plt.figure() 
    ax = plt.axes(xlim=(ptc.xax[0], ptc.xax[-1]), ylim=(-1.5, 1.5)) 
    liner, = ax.plot([], [], lw=0.5) 
    linei, = ax.plot([], [], lw=0.5) 
    linep, = ax.plot([], [], lw=0.7)    
    linev, = ax.plot([], [], lw=2)

    def init(): 
        liner.set_data([], []) 
        linei.set_data([], []) 
        linep.set_data([], [])

        linev.set_data([], [])

        return liner, linei, linep, linev
 
    def animate(i):
        global ptc
        ptc.update()

        stater = np.real(ptc.psi_x) 
        statei = np.imag(ptc.psi_x) 
        p = stater**2 + statei**2 
        liner.set_data(ptc.xax, stater) 
        linei.set_data(ptc.xax, statei) 
        linep.set_data(ptc.xax, p)

        linev.set_data(ptc.xax, ptc.poten)
         
        return liner, linei, linep, linev
 
    t0 = time.time() 
    animate(0) 
    t1 = time.time() 
    intv = int(1000 * ptc.dt - (t1-t0)) 
    frms = (tot_time*1000)//intv 
    print(intv) 
 
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=frms, interval=intv, blit=True) 
    plt.show() 


def fourier(x, psi_x):
    assert len(x) == len(psi_x)
    N = len(x)
    dx = x[1] - x[0]
    k0 = -np.pi/dx
    sr2pi = np.sqrt(2*np.pi)
    dk = 2*np.pi/(N*dx)

    k = np.arange(0, N) * dk + k0
    psi_k_1 = fft((dx/sr2pi)*psi_x*np.exp(-1j*k0*x))
    psi_k = psi_k_1*np.exp(1j*(k-k0)*x[0])
    return k, psi_k,

def ifourier(k, psi_k, x):
    assert len(k) == len(psi_k) and len(k) == len(x)
    sr2pi = np.sqrt(2*np.pi)
    dx = x[1]-x[0]
    psi_x_1 = ifft(psi_k*np.exp(-1j*(k-k[0]*x)))
    psi_x = psi_x_1 * sr2pi / dx * np.exp(1j*k[0]*x)
    return psi_x

class ParticleInABox:
    def __init__(self, x, psi_x0, m, pot=0, dt=0.01):
        self.xax = x
        self.psi_x = psi_x0
        self.N = len(self.xax)

        self.kax, tmp = fourier(self.xax, self.psi_x)
        
        assert len(x) == len(psi_x0) and self.N >= 512

        self.mass = m
        self.poten = wall(self.N*3//100, self.xax)
        self.poten += pot
        self.dt = dt
        self.time = 0
        
        self.ux_half = np.exp(-1j*self.poten*self.dt/(2*hbar))
        self.ux = self.ux_half * self.ux_half
        self.uk = np.exp(-1j*self.kax*self.kax*hbar*self.dt/(2*self.mass))

        
    def update(self):
        self.psi_x = self.psi_x * self.ux_half
        ktmp, psi_k = fourier(self.xax, self.psi_x)
        psi_k = psi_k * self.uk
        self.psi_x = ifourier(self.kax, psi_k, self.xax)
        self.psi_x = self.psi_x * self.ux_half
    
        self.time += self.dt
        
x = np.linspace(-100, 100, 2**12)
psi_x = gauss_packet(x, 0, delta=2)

ptc = ParticleInABox(x, psi_x, 10000)

animator(ptc, 10)




