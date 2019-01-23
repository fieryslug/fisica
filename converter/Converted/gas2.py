from vpython import *
from time import clock
#from vpython.graph import *
from random import random

# A model of an ideal gas with hard-sphere collisions
# Program uses Numeric Python arrays for high speed computations

win=500

# change this to have more or fewer atoms
Natoms = 50  

# Typical values
# container is a cube L on a side
L = 1. 
# color of edges of container
gray = vector(0.7,0.7,0.7) 
# radius of lines drawn on edges of cube
Raxes = 0.005 
# helium mass
Matom = 4E-3/6E23 
# wildly exaggerated size of helium atom
Ratom = 0.03 
# Boltzmann constant
k = 1.4E-23 
# around room temperature
T = 300. 
dt = 1E-5

scene = canvas(title="Gas", width=win, height=win, x=0, y=0, range=L, center=vector(L/2.,L/2.,L/2.))

# binning for v histogram
deltav = 100. 
vdist = graph(x=0, y=win, ymax = Natoms*deltav/1000., width=win, height=win/2, xtitle='v', ytitle='dN')
theory = gcurve(color=color.cyan)
observation = ghistogram(bins=arange(0.,3000.,deltav), accumulate=1, average=1, color=color.red)

dv = 10.
# theoretical prediction
for v in arange(0.,3001.+dv,dv): 
    theory.plot(pos=(v, (deltav/dv)*Natoms*4.*pi*((Matom/(2.*pi*k*T))**1.5) *exp((-0.5*Matom*v**2)/(k*T))*v**2*dv))

xaxis = curve(pos=[(0,0,0), (L,0,0)], color=gray, radius=Raxes)
yaxis = curve(pos=[(0,0,0), (0,L,0)], color=gray, radius=Raxes)
zaxis = curve(pos=[(0,0,0), (0,0,L)], color=gray, radius=Raxes)
xaxis2 = curve(pos=[(L,L,L), (0,L,L), (0,0,L), (L,0,L)], color=gray, radius=Raxes)
yaxis2 = curve(pos=[(L,L,L), (L,0,L), (L,0,0), (L,L,0)], color=gray, radius=Raxes)
zaxis2 = curve(pos=[(L,L,L), (L,L,0), (0,L,0), (0,L,L)], color=gray, radius=Raxes)

Atoms = []
colors = [color.red, color.green, color.blue,
          color.yellow, color.cyan, color.magenta]
poslist = []
plist = []
mlist = []
rlist = []

for i in range(Natoms):
    Lmin = 1.1*Ratom
    Lmax = L-Lmin
    x = Lmin+(Lmax-Lmin)*random()
    y = Lmin+(Lmax-Lmin)*random()
    z = Lmin+(Lmax-Lmin)*random()
    r = Ratom
    Atoms = Atoms+[sphere(pos=vector(x,y,z), radius=r, color=colors[i % 6])]
    mass = Matom*r**3/Ratom**3
# average kinetic energy p**2/(2mass) = (3/2)kT
    pavg = sqrt(2.*mass*1.5*k*T) 
    theta = pi*random()
    phi = 2*pi*random()
    px = pavg*sin(theta)*cos(phi)
    py = pavg*sin(theta)*sin(phi)
    pz = pavg*cos(theta)
    poslist.append((x,y,z))
    plist.append((px,py,pz))
    mlist.append(mass)
    rlist.append(r)

pos = array(poslist)
p = array(plist)
m = array(mlist)
# Numeric Python: (1 by Natoms) vs. (Natoms by 1)
m.shape = vector(Natoms,1,0) 
radius = array(rlist)

t = 0.0
Nsteps = 0
# initial half-step
pos = pos+(p/m)*(dt/2.) 
time = clock()

while 1:
    observation.plot(data=mag(p/m))

    # Update all positions
    pos = pos+(p/m)*dt

# all pairs of atom-to-atom vectors
    r = pos-pos[:,NewAxis] 
# atom-to-atom scalar distances
    rmag = sqrt(add.reduce(r*r,-1)) 
    hit = less_equal(rmag,radius+radius[:,NewAxis])-identity(Natoms)
# i,j encoded as i*Natoms+j
    hitlist = sort(nonzero(hit.flat)).tolist() 

    # If any collisions took place:
    for ij in hitlist:
# decode atom pair
        i, j = divmod(ij,Natoms) 
# remove symmetric j,i pair from list
        hitlist.remove(j*Natoms+i) 
        ptot = p[i]+p[j]
        mi = m[i,0]
        mj = m[j,0]
        vi = p[i]/mi
        vj = p[j]/mj
        ri = Atoms[i].radius
        rj = Atoms[j].radius
        a = mag(vj-vi)**2
# exactly same velocities
        if a == 0: continue 
        b = 2*dot(pos[i]-pos[j],vj-vi)
        c = mag(pos[i]-pos[j])**2-(ri+rj)**2
        d = b**2-4.*a*c
# something wrong; ignore this rare case
        if d < 0: continue 
# t-deltat is when they made contact
        deltat = (-b+sqrt(d))/(2.*a) 
# back up to contact configuration
        pos[i] = pos[i]-(p[i]/mi)*deltat 
        pos[j] = pos[j]-(p[j]/mj)*deltat
        mtot = mi+mj
# transform momenta to cm frame
        pcmi = p[i]-ptot*mi/mtot 
        pcmj = p[j]-ptot*mj/mtot
        rrel = norm(pos[j]-pos[i])
# bounce in cm frame
        pcmi = pcmi-2*dot(pcmi,rrel)*rrel 
        pcmj = pcmj-2*dot(pcmj,rrel)*rrel
# transform momenta back to lab frame
        p[i] = pcmi+ptot*mi/mtot 
        p[j] = pcmj+ptot*mj/mtot
# move forward deltat in time
        pos[i] = pos[i]+(p[i]/mi)*deltat 
        pos[j] = pos[j]+(p[j]/mj)*deltat
 
    # Bounce off walls
# walls closest to origin
    outside = less_equal(pos,Ratom) 
    p1 = p*outside
# force p component inward
    p = p-p1+abs(p1) 
# walls farther from origin
    outside = greater_equal(pos,L-Ratom) 
    p1 = p*outside
# force p component inward
    p = p-p1-abs(p1) 

    # Update positions of display objects
    for i in range(Natoms):
        Atoms[i].pos = pos[i]

    Nsteps = Nsteps+1
    t = t+dt

    if Nsteps == 50:
        print('%3.1f seconds for %d steps with %d Atoms' % (clock()-time, Nsteps, Natoms))
##    rate(30)

