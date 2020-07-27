###
# TODO List
###
# ☒ add diploid stuff
# ☐ fix animation
#       * something's broken with the animation, and I'm not really sure what
#       * I think it has to do with imgs.set_array
# ☐ add a field in the chromosome saying which index of 
#   morphogen it codes for.
#       * Right now, a haploid genome always codes for morphogens 0 to N-1
#         they can't skip morphogen 2, for example
# ☐ add some way to linearize each chromosome
# ☐ support reproduction between two diploid genomes, including crossover
#   and chromosome duplication / deletion
# ☐ add cell differentiation
#       * Based on the levels of particular morphogens, each x, y location will be assigned a cell type
#       * This includes:
#            + muscle
#            + bone
#            + filler tissue
#            + not a cell (counts as empty space)
# ☐ add motion simulation
#       * A very low-tech fluid simulation
#       * assumes water molecules don't move in space
#       * applies force to creature based on how many gridpoints it hits between shape changes and what direction those points were hit from
#       * changes creature's shape by contracting muscle cells
#            + I'm not totally sure how to compute the creature's new shape yet, but I'm thinking of using something like photoshop's distort tool
#            + I also have no idea how or when the creature's muscles will activate either

#from benmaier
from benmaier_reactionDiffusion.gray_scott_static import get_initial_A_and_B
from benmaier_reactionDiffusion.gray_scott_static import update
import matplotlib.pyplot as pl
import numpy as np
import random

# animation imports from benmaier
import matplotlib.animation as animation
from matplotlib.colors import Normalize


# modified from benmaier for drawing multiple pairs
def draw(pairs):
    """return the matplotlib artists for animation"""
    imgs = []
    
    i = 0
    # create one row per pair, where each row has 2 columns
    fig, ax = pl.subplots(len(pairs),2,figsize=(5.65,3), squeeze=False)
    for A, B in pairs:
        imA = ax[i][0].imshow(A, animated=True,vmin=0,cmap='Greys')
        imB = ax[i][1].imshow(B, animated=True,vmax=1,cmap='Greys')
        imgs.append([imA, imB])
        
        ax[i][0].axis('off')
        ax[i][1].axis('off')
        ax[i][0].set_title('A'+str(i))
        ax[i][1].set_title('B'+str(i))
        i += 1

    return fig, imgs


# modified from benmaier
def hardcodedMain():
    # =========== define model parameters ==========

    # update in time
    delta_t = 1.0

    # Diffusion coefficients (diffusion speeds)
    DA = 0.16
    DB = 0.08

    # define create/decay rates
    f = 0.060
    k = 0.062

    # grid size
    N = 200

    # intialize the chemical concentrations
    A, B = get_initial_A_and_B(N)
    C, D = get_initial_A_and_B(N)
    E, F = get_initial_A_and_B(N)

    N_simulation_steps = 1000#10000
    for step in range(N_simulation_steps):
        A, B = update(A, B, DA, DB, f, k, delta_t)
        C, D = update(C, D, DA, DB, A*0.1, k, delta_t)
        E, F = update(E, F, DA, DB, A*A, k, delta_t)
    print(A)

    draw( [(A, B), (C, D), (E, F)] )

    # show the result
    pl.show()
    
  
  
def defaultChromosome():
    return                        \
    {                             \
        "f": {"-":0.060},         \
        "k": {"-":0.062},         \
        "DA":{"-":0.16},          \
        "DB":{"-":0.08}           \
    }                             \

def defaultDependancyChromosome():
    return                        \
    {                             \
        "f":                      \
            {                     \
                "-":0.0,          \
                "A0":0.1          \
            },                    \
        "k": {"-":0.060},         \
        "DA":{"-":0.16},          \
        "DB":{"-":0.08}           \
    }                             \


def randomChromosome(N):
    chromosome = defaultChromosome()
    
    for geneName in chromosome:
        gene = chromosome[geneName]
        
        # 50%   - this gene has no dependancies
        # 25%   - this gene has 1 dependancy
        # 12.5% - this gene has 2 depenancies ...
        for i in range(N-1):
            if random.random() < 0.5:
                break
            # to reduce headaches, any gene with a dependancy automatically has only that dependancy
            gene["-"] = 0 
            
            # 80% - this dependancy is on A
            # 20% - this dependancy is on B
            dependancy  = "A" if random.random() < 0.8 else "B"
            dependancy += str(random.randint(0, N-1))
            amt = random.random()*0.1
            
            gene[dependancy] = amt
    return chromosome


def randomGenome(N):
    genome = [defaultChromosome()]
    for i in range(N-1):
        genome.append(randomChromosome(N))
    return genome


# if a key is "-" then it codes for just a scalar, not a dependancy coefficient
def decode(gene, morphogens): 
    result = 0
    for key in gene:
        if key == '-': # not a dependancy
            result += gene[key]
        else:
            index = int(key[1:])
            result += gene[key]*morphogens[index][0 if key[0]=="A" else 1]
            
    return result


def __simulationStep(chromosomes, morphogens, delta_t, gridSize):
    newMorphogens = []
    # calculate parameters (important for dependancies)
    fValues = []
    kValues = []
    DAValues = []
    DBValues = []
    for chromosome in chromosomes:
        fValues.append(decode(chromosome["f"], morphogens))
        kValues.append(decode(chromosome["k"], morphogens))
        DAValues.append(decode(chromosome["DA"], morphogens))
        DBValues.append(decode(chromosome["DB"], morphogens))
    
    
    for morphogen, f, k, DA, DB in zip(morphogens, fValues, kValues, DAValues, DBValues):
        newMorphogens.append(tuple(update(*morphogen, DA, DB, f, k, delta_t)))
    
    # in case the given set of chromosomes don't cover every morphogen
    for i in range(len(morphogens) - len(newMorphogens)):
        newMorphogens.append(np.zeros((gridSize,gridSize)))
    
    return newMorphogens
    
    
def simulate(chromosomesF, chromosomesM, gridSize=200, N_simulation_steps=1000, delta_t=1.0):
    numChromosomes = max(len(chromosomesF), len(chromosomesM))
    morphogens = [np.asarray(get_initial_A_and_B(gridSize)) for i in range(numChromosomes)]
    
    morphogens = __computeSimulation(morphogens, chromosomesF, chromosomesM, delta_t, gridSize, N_simulation_steps)
    
    draw( morphogens )

    # show the result
    pl.show()
    return morphogens


def __computeSimulation(morphogens, chromosomesF, chromosomesM, delta_t, gridSize, N_simulation_steps):
    for step in range(N_simulation_steps):
        morphogensF = __simulationStep(chromosomesF, morphogens, delta_t, gridSize)
        morphogensM = __simulationStep(chromosomesM, morphogens, delta_t, gridSize)
        
        morphogens = [0.5*np.asarray(f) + 0.5*np.asarray(m) for f, m in zip(morphogensF, morphogensM)]
    return morphogens
    
 
# modified from benmaier for multiple pairs
def updatefig(frame_id, *args):
    """Takes care of the matplotlib-artist update in the animation"""
    #print(*args)
    
    imgs = args[0]
    simulationArgs = args[1:]
    morphogens = __computeSimulation(*simulationArgs)
    #args[1] = morphogens
    for i in range(len(args[1])):
        args[1][i] = morphogens[i]

    #print(str(frame_id) + " ==========*************=============")

    for i in range(len(morphogens)):
        A = morphogens[i][0]
        B = morphogens[i][1]
        #print(A.shape)
        #print(A)
        #print("-==--=-=-=-=-=-=-=-=-=-")
        
        # update the frame
        imgs[i][0].set_array(A)
        imgs[i][1].set_array(B)
        
        # renormalize the colors
        imgs[i][0].set_norm(Normalize(vmin=np.amin(A),vmax=np.amax(A)))
        imgs[i][1].set_norm(Normalize(vmin=np.amin(B),vmax=np.amax(B)))
        
    # return the updated matplotlib objects
    return (i for pair in imgs for i in pair)


# modified from benmaier 
def animate(chromosomesF, chromosomesM, gridSize=200, updates_per_frame=10, delta_t=1.0):
    numChromosomes = max(len(chromosomesF), len(chromosomesM))
    morphogens = [np.asarray(get_initial_A_and_B(gridSize)) for i in range(numChromosomes)]
    
    
    fig, imgs = draw(morphogens)
    # these are the arguments which have to passed to the update function
    animation_arguments = (imgs, morphogens, chromosomesF, chromosomesM, delta_t, gridSize, updates_per_frame)

    # start the animation
    ani = animation.FuncAnimation(fig, #matplotlib figure
                                  updatefig, # function that takes care of the update
                                  fargs=animation_arguments, # arguments to pass to this function
                                  interval=1, # update every `interval` milliseconds
                                  blit=True, # optimize the drawing update 
                                 )

    # show the animation
    pl.show()


if __name__=="__main__":
    #genome =                \
    #[                       \
        #defaultGeneSet(),   \
        #defaultDependancyGeneSet() \
    #]
    #simulate(genome)
    
    
    import pprint
    pp = pprint.PrettyPrinter(indent=4)
    
    chromosomesF = randomGenome(4)
    chromosomesM = randomGenome(4)
    
    pp.pprint(chromosomesF)
    pp.pprint(chromosomesM)
    
    #simulate(chromosomesF, chromosomesM)
    animate(chromosomesF, chromosomesM)
