import math
from matplotlib import pyplot

h_ = 1
pi = math.pi
n = 2
a = 6

def prob_p(p):
    return a*pi*(n**2)*(2-2*math.cos(p*a/h_ - n*pi)) / (h_ * ((n*pi)**2 - (a*p/h_)**2)**2)



def calc_ppave():
    p0 = -100
    pmax = 100
    inc = 0.001

    pp_ave = 0
    
    while p0 < pmax:
        pp_ave += prob_p(p0)*inc
        p0 += inc
        print(p0)

    print(pp_ave)

def plot_prob():
    p0 = -10
    pmax = 10

    pplot = []
    probplot = []

    inc = 0.01

    while p0 <= pmax:

        pplot.append(p0)
        probplot.append(prob_p(p0))
        p0 += inc


    pyplot.plot(pplot, probplot)
    pyplot.show()

calc_ppave()
