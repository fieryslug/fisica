from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np
import common_1 as mq


fig = plt.figure()
ax = plt.axes(xlim=(-8, 8), ylim=(-3, 3))
liner, = ax.plot([], [], lw=0.5)
linei, = ax.plot([], [], lw=0.5)
linep, = ax.plot([], [], lw=0.7)

def main():
    f = 1
    x = np.linspace(-8, 8, 1000)
    ptc = mq.Particle(1, -8, 8)
    ptc.poten = -f * x
    ptc.state = mq.gauss_packet(x)

    mq.animator(ptc, 0.01, 10)
        
    



if __name__ == '__main__':
    main()
