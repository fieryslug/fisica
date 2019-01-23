

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
    pygame.display.set_caption("Griffiths")
    screen.fill((0, 0, 0))

    particles = []
    
    running = True

    currtime = time.time()

    N = len(particles)

    while running:
        if time.time() - currtime >= dt:
            currtime = time.time()
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    running = False
            for ptc in particles:
                ptc.update()
            screen.fill(background)
             

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



        
