import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import plotter as pltr
from scipy.fftpack import fft, ifft, fftshift
import time

hbar = 1.0

def gauss_packet(x, x0=0, p0=0, delta=1):
    return np.exp(1j*p0*(x-x0)/hbar) * np.exp(-(x-x0)**2/(2*delta**2)) / ((np.pi*delta**2)**(1/4))

def plane_wave(x, p0):
    return np.exp(1j*p0*x/hbar)

def wall(width, x):
    y = 0 * x
    N = len(y)
    for i in range(width):
        y[i] = 1e6
        y[-1-i] = 1e6
    return y

def theta(x):
    return np.heaviside(x, 0)

def delta(x, x0=0):
    sig = 0.05
    return 1 / (sig * np.sqrt(2*np.pi)) * np.exp(-0.5 * ((x-x0)/sig)**2)
    



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
        self.dx = self.xax[1] - self.xax[0]
        self.psi_x = psi_x0
        self.N = len(self.xax)

        self.kax, tmp = fourier(self.xax, self.psi_x)
        self.kax = np.linspace(-1, 1, self.N) * (np.pi/self.dx)
        self.dk = self.kax[1] - self.kax[0]

        assert len(x) == len(psi_x0) and self.N >= 512

        self.mass = m
        self.poten = wall(self.xax.size//30, self.xax)
        #self.poten = 0
        self.poten = np.maximum(self.poten, pot)
        self.dt = dt
        self.time = 0
        
        #self.ux_half = np.exp(-1j*self.poten*self.dt/(2*hbar))
        #self.ux = self.ux_half * self.ux_half
        #self.uk = np.exp(-1j* self.kax * self.kax * hbar * self.dt/ (2*self.mass))

    def xtok(self):
        psi_x_mod = self.dx / np.sqrt(np.pi*2) *np.exp(1j*self.kax[0]*self.xax)*self.psi_x
        psi_k_mod = fft(psi_x_mod)
        psi_k = psi_k_mod * np.exp(1j*(self.kax-self.kax[0])*self.xax[0])
        return psi_k

    def ktox(self, psi_k):
        psi_k_mod = psi_k * np.exp(-1j*(self.kax-self.kax[0])*self.xax[0])
        psi_x_mod = ifft(psi_k_mod)
        self.psi_x = psi_x_mod * np.sqrt(2*np.pi) / self.dx * np.exp(-1j*self.kax[0]*self.xax)

    @property
    def prob(self):
        return np.abs(self.psi_x) ** 2
        
    def update(self):
        
        self.psi_x *= np.exp(-0.5j/hbar * self.poten * self.dt)
        psi_k = self.xtok()
        psi_k *= np.exp(-1j*self.kax*self.kax*hbar/(2*self.mass) * self.dt)
        self.ktox(psi_k)
        self.psi_x *= np.exp(-0.5j/hbar * self.poten * self.dt)
        self.time += self.dt

def main():
    
   
    pltr.cplot(ptc.kax, ptc.xtok())
    fig = plt.figure()
    ax = plt.axes(xlim = (ptc.xax[0], ptc.xax[-1]), ylim = (-2, 2))
    liner, = ax.plot([], [], lw=1)
    linei, = ax.plot([], [], lw=1)
    linep, = ax.plot([], [], lw=1)
    linev, = ax.plot([], [], lw=1)

    def init(): 
        liner.set_data([], []) 
        linei.set_data([], []) 
        linep.set_data([], [])
        linev.set_data([], [])

        return (liner, linei, linep, linev)
     
    def animate(i):
        global ptc
        ptc.update()

        stater = np.real(ptc.psi_x) 
        statei = np.imag(ptc.psi_x) 
        p = stater**2 + statei**2 
        liner.set_data(ptc.xax, stater) 
        linei.set_data(ptc.xax, statei) 
        linep.set_data(ptc.xax, p)

        linev.set_data(ptc.xax, ptc.poten/50)
        print(str(np.around(ptc.time, 2)))
             
        return (linep, linev)

    t0 = time.time()
    animate(0)
    t1 = time.time()
    intv = 1000 * ptc.dt - (t1-t0)

    anim = animation.FuncAnimation(fig, animate, init_func=init, interval=intv, blit=True, frames=300)



if __name__ == '__main__':
    main()


