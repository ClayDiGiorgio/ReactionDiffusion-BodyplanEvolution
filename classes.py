import random

class Constants:
    # below should probably be on (0, 1]
    diffusionRateP = 1
    diffusionRateS = 0.5 # should be lower than diffusionRateP
    pAddRate = 0.055   # these two are usually what changes to affect the pattern
    sDecayRate = 0.062 #
    sAddRate = 0.01
    
    reproduceLevel = 0.5
    dieLevel = 0.5
    
    neighborDistance = 10
    
    cellLimit = 100
    
    friction = 0.1
    

class Utilities:
    def unitizeVector(vector):
        mag = (vector[0]**2 + vector[1]**2)**0.5
        if mag == 0:
            return (0, 0)
        return (vector[0]/mag, vector[1]/mag)

    def setVectorMagnitude(vector, mag):
        vec = Utilities.unitizeVector(vector)
        return (mag*vec[0], mag*vec[1])
    
    def vectorDistance(v, w):
        return ((v[0]-w[0])**2 + (v[1]-w[1])**2)**0.5
    
    def vectorScale(v, s):
        return (v[0]*s, v[1]*s)
    
    def vectorAdd(v, w):
        return (v[0]+w[0], v[1]+w[1])
     
    def vectorSubtract(v, w):
        return (v[0]-w[0], v[1]-w[1])
    
    def randomDirectionVector():
        return (random.uniform(-1, 1), random.uniform(-1, 1))
     

class Organism:
    def __init__(self):
        stemCell = Cell(self, (0,0))
        stemCell.p = 1
        
        self.cells = [stemCell]
        self.step = 0
    
    
    def update(self, deltaT):
        for c in self.cells:
            c.reactToMorphogenLevels()
            c.updateMorphogenLevels(deltaT)
            c.updateLocation(deltaT)
            #NOTE: these should probably be done in separate loops
            #NOTE: should probably update all cells' morphogen levels based on current timestep and then commit all changes together after calculations are done
        self.step += 1
    
    def spawnCell(self, parent, direction):
        print("step "+str(self.step)+": spawn!")
        if len(self.cells) > Constants.cellLimit:
            return
        
        if direction == (0,0):
            #direction = (random.uniform(-1,1), random.uniform(-1,1))
            direction = (0, -1)
        offset = Utilities.setVectorMagnitude(direction, Constants.neighborDistance)
        location = (parent.location[0]+offset[0], parent.location[1]+offset[1])
        
        cell = Cell(self, location)
        self.cells.append(cell)
        #parent.neighbors.append(cell)
        
        for c in self.cells:
            if Utilities.vectorDistance(c.location, location) <= Constants.neighborDistance:
                c.neighbors.append(cell)
                cell.neighbors.append(c)
        
    
    def removeCell(self, cell):
        for c in cell.neighbors:
            try:
                cell.neighbors.remove(cell)
            except:
                anErrorHappenedThatIDontCareAbout=True
        #cells.remove(cell)
    
    
    

class Cell:
    #organism
    
    #x, y
    #list connected cells
    #p level
    #s level
    
    #p source direction
        #| caluclated from the average relative location of all neighbors, scaled by their relative levels of p
    #s source direction
        #| caluclated from the average relative location of all neighbors, scaled by their relative levels of s
    
    #*note: p flow direction would be -(p source direction)
        
    #eg: 
        #(where P = (neighbor.p-self.p) =  1 
         #and   p = (neighbor.p-self.p) = -1 )
        
        #P p 0
        #0 O 0
        #P 0 0
        
        #p direction would be AVG(<-1, -1> + <-1, 1> - <0, 1>) = <-0.66, -0.33>
    
    
    def __init__(self, organism, location):
        #assert(type(organism) == Organism)
        #assert(type(location) == tuple)
        #assert(len(location) == 2)
        #assert(type(location[0]) == type(0.0))
        #assert(type(location[1]) == type(0.0))
        
        self.organism = organism
        self.location = location
        self.neighbors = []
        self.p = 0
        self.s = 0
        self.pSourceDirection = (0, 0)
        self.sSourceDirection = (0, 0)
        
        self.velocity = (0, 0)
        
    
    def updateLocation(self, deltaT):
        force = (0, 0)
        
        if self.velocity != (0, 0):
            friction = Utilities.setVectorMagnitude(self.velocity, -Constants.friction)
            force = Utilities.vectorAdd(force, friction)
        
        for c in self.neighbors:
            if c == self:
                continue
            distance = Utilities.vectorDistance(self.location, c.location)
            if distance > Constants.neighborDistance:
                mag = Cell.connectionForce(distance)
                Dir = Utilities.vectorSubtract(c.location, self.location)
                if Dir == (0, 0):
                    Dir = Utilities.randomDirectionVector()
                connectForce = Utilities.setVectorMagnitude(Dir, mag)
                force = Utilities.vectorAdd(force, connectForce)
        
        for c in self.organism.cells:
            if c == self:
                continue
            distance = Utilities.vectorDistance(self.location, c.location)
            if distance < Constants.neighborDistance:
                mag = Cell.pushingForce(distance)
                Dir = Utilities.vectorSubtract(c.location, self.location)
                if Dir == (0, 0):
                    Dir = Utilities.randomDirectionVector()
                pushForce = Utilities.setVectorMagnitude(Dir, mag)
                force = Utilities.vectorAdd(force, pushForce)
        
        
        self.velocity = Utilities.vectorAdd(self.velocity, force)
        self.location = Utilities.vectorAdd(self.location, Utilities.vectorScale(self.velocity, deltaT))
        
        #for c in self.neighbors:
            #if Utilities.vectorDistance(self.location, c.location) > Constants.neighborDistance:
                #self.neighbors.remove(c)
                #c.neighbors.remove(self)

    
    def connectionForce(distance):
        distance = distance-Constants.neighborDistance
        
        scaling = 10
        maxForce = 10
        return maxForce - 1/(scaling*distance + 1)
    
    def pushingForce(distance):
        distance = distance-Constants.neighborDistance
        
        scaling = 10
        return 0.1 - 1/(scaling*distance +10)
    
    
    # formulas modified from https://www.karlsims.com/rd.html
    def updateMorphogenLevels(self, deltaT):
        neighbors = self.neighbors
        p = self.p
        s = self.s
        
        averageNeighborP = 0
        averageNeighborS = 0
        if len(neighbors) > 0:
            averageNeighborP = sum(c.p for c in neighbors)/len(neighbors)
            averageNeighborS = sum(c.s for c in neighbors)/len(neighbors)
        
        self.p =                                                   \
            self.p                                                 \
            + (                                                    \
                Constants.diffusionRateP * (averageNeighborP-p)    \
                - p*s*s                                            \
                + Constants.pAddRate*(1-p)                         \
            )*deltaT                                               
            
        self.s =                                                   \
            self.s                                                 \
            + (                                                    \
                Constants.diffusionRateS * (averageNeighborS-s)    \
                + p*s*s                                            \
                - (Constants.sDecayRate+Constants.pAddRate)*s      \
                + Constants.sAddRate*(1-s)
            )*deltaT
         
        pSourceDirection =                                         \
        (                                                          \
            sum((c.p-self.p)*(c.location[0]-self.location[0]) for c in neighbors),     \
            sum((c.p-self.p)*(c.location[1]-self.location[1]) for c in neighbors)      \
        )
        
        sSourceDirection =                                         \
        (                                                          \
            sum((c.s-self.s)*(c.location[0]-self.location[0]) for c in neighbors),     \
            sum((c.s-self.s)*(c.location[1]-self.location[1]) for c in neighbors)      \
        )
        
            
    def reactToMorphogenLevels(self):
        if self.p > Constants.reproduceLevel:
            #spawn new cell in p flow direction
            #add cell to neighbor list of all cells within Constants.neighborDistance distance of it
                #and also add those cells to its neighbor distance
            self.organism.spawnCell(self, Utilities.vectorScale(self.pSourceDirection, -1))    
            
        if self.s > Constants.dieLevel:
            #delete this cell and remove it from the neighbor list of all its neighbors
            self.organism.removeCell(self)
            
    #updateLocation()
        #maintain a constant distance from all neighbors (defined by Constants.neighborDistance)
        #repel from 
        
    def __str__(self):
        return "Cell: p="+str(self.p)+", s="+str(self.s)
