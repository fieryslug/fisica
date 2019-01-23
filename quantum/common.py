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
        y[i] = 1e6
        y[N-1-i] = 1e6
    return y

def animator(dt, tot_time, x, psi_x): 
    fig = plt.figure() 
    ax = plt.axes(xlim=(-100, 100), ylim=(-1.5, 1.5)) 
    liner, = ax.plot([], [], lw=0.5) 
    linei, = ax.plot([], [], lw=0.5) 
    linep, = ax.plot([], [], lw=0.7)    


    def init(): 
        liner.set_data([], []) 
        linei.set_data([], []) 
        linep.set_data([], []) 
        return liner, linei, linep, 
 
    def animate(i):
        time.sleep(1)
        global x, psi_x
        m = 1
        k, psi_k = fourier(x, psi_x)
        psi_k += dt * (k*k*hbar/(2*m*1j)) * psi_k
        psi_x = ifourier(k, psi_k, x)

        stater = np.real(psi_x) 
        statei = np.imag(psi_x) 
        p = stater**2 + statei**2 
        liner.set_data(x, stater) 
        linei.set_data(x, statei) 
        linep.set_data(x, p) 
         
        return liner, 
 
    t0 = time.time() 
    animate(0) 
    t1 = time.time() 
    intv = int(1000 * dt - (t1-t0)) 
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

x = np.linspace(-100, 100, 2**10)
psi_x = gauss_packet(x, 5)
animator(0.01, 10, x, psi_x)




