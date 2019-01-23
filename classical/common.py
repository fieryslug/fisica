import pygame

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


