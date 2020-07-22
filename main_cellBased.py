from classes import Organism
from classes import Constants
import pygame
import time

organism = Organism()

width = 500
height = 500

window = pygame.display.set_mode((width, height))
organismCellView = window

timeStep = 1
cellRadius = int(Constants.neighborDistance/2.0)

done = False
while not done:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            pygame.quit()
            done = True
            break
        
    
    organism.update(timeStep)
    
    
    organismCellView.fill((100,200,150,))
    for c in organism.cells:
        center = (int(c.location[0]+width/2), int(c.location[1]+height/2))
        color = (int(255*c.p),)*3
        pygame.draw.circle(organismCellView, color, center, cellRadius)
        
        
    pygame.display.flip()

    time.sleep(0.1)
