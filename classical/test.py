import pygame
import math
import time
import hamilton

pixel_per_meter = 10
width = 800
height = 600
v_inf = 1e4
width_m = width/pixel_per_meter
hegiht_m = height/pixel_per_meter
mass0 = 1


def main():
    dims = (width, height)
    background = (32, 32, 32)
    screen = pygame.display.set_mode(dims)
    pygame.display.set_caption("h")
    screen.fill((0, 0, 0))

    particles = []
    
    P = [10, 10]
    Q = [10, 10]

    i = 0
    while i < len(P):
        particle0 = Particle(screen, Vec2(Q[i], Q[i+1]), 5)
        particle0.p.x = P[i]
        particle0.p.y = P[i+1]
        particles.append(particle0)
        i += 2

    running = True

    currtime = time.time()
    dt = 0.001
    dP = 0.0001
    dQ = 0.0001

    N = len(particles)

    while running:
        if time.time() - currtime >= dt:
            currtime = time.time()
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    running = False
            screen.fill(background)
                    
            (P, Q) = hamilton.update(P, Q, Hamilton, dt)
            copyPhase(P, Q, particles)
            print(P[0] ** 2 + P[1] ** 2)
            for i in range(N):
                #particles[i].update(dt)
                particles[i].display()
                pygame.display.flip()

def copyPhase(P, Q, particles):
    i = 0
    j = 0
    while i < len(P):
        particles[j].pos.x = Q[i]
        particles[j].pos.y = Q[i+1]
        particles[j].p.x = P[i]
        particles[j].p.y = P[i+1]
        i += 2
        j += 1

def Hamilton(P, Q):
    T = 0
    for pi in P:
        T += (pi**2) / (2*mass0)
    V = 0
    for i in range(len(Q)):
        if i % 2 == 0:#qx
            V += v_inf*(unit_step(-Q[i]) + unit_step(Q[i]-width/10))

        else:#qy
            V += v_inf*(unit_step(-Q[i]) + unit_step(Q[i]-height/10))

    return T + V

def unit_step(x):
    return (1 + math.tanh(0.5*x)) / 2

class Particle:
    def __init__(self, screen, pos, rad):
        self.screen = screen
        self.pos = pos
        self.p = Vec2(0, 0)
        self.radius = rad
        self.mass = mass0

    @property
    def vel(self):
        return Vec2(self.p.x/self.mass, self.p.y/self.mass)

    def display(self):
        pygame.draw.circle(self.screen, (192, 192, 192), (int(self.pos.x*pixel_per_meter), int(self.pos.y*pixel_per_meter)), self.radius, 2)

    def update(self, dt):
        self.pos.x += self.vel.x * dt
        self.pos.y += self.vel.y * dt


class Vec2:
    def __init__(self, x, y):
        self.x = x
        self.y = y

if __name__ == '__main__':
    main()
