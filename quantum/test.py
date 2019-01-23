import common_1 as mq
import numpy as np
import plotter as pltr

ptc = mq.Particle(1, -5, 5, pres=1000)
ptc.state = mq.gauss_packet(ptc.x)

f = 4 * ptc.x
g = ptc.ddx(ptc.ddx(f))

for i in range(4):
    ptc.update(0.01)

#pltr.cplot(ptc.x, ptc.ddx(ptc.state))
pltr.cplot(ptc.x, ptc.Hamilton())

