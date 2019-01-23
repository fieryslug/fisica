import hamilton
from matplotlib import pyplot

def main():
    tottime = 10
    dt = 0.0001
    P = [0]
    Q = [2]
    
    Pplot = []
    Qplot = []
    tplot = []
    t = 0
    while t <= tottime:
        (P, Q) = hamilton.update(P, Q, Hamilton, dt)
        t += dt
        Pplot.append(P[0])
        Qplot.append(Q[0])
        tplot.append(t)

    #pyplot.plot(tplot, Pplot, 'r', tplot, Qplot, 'g')
    pyplot.plot(Pplot, Qplot)
    pyplot.show()
        

def Hamilton(P, Q):
    k = 1
    m = 1
    p = P[0]
    q = Q[0]
    return p**2/(2*m) + (1/2)*k*q**2



if __name__ == '__main__':
    main()
