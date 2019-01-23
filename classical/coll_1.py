import pygame
import math
import time
import hamilton
import common

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
    

    running = True
    currtime = time.time()


    dt = 0.001
    N = len(particles)

    while running:
        if time.time() - currtime >= dt:
            currtime = time.time()
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    running = False
            screen.fill(background)


