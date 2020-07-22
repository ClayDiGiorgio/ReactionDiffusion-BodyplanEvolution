#from benmaier_reactionDiffusion.gray_scott_static import draw
from benmaier_reactionDiffusion.gray_scott_static import get_initial_A_and_B
from benmaier_reactionDiffusion.gray_scott_static import update
import matplotlib.pyplot as pl
import numpy as np


# modified from benmaier for drawing multiple pairs
def draw(pairs):
    """return the matplotlib artists for animation"""
    i = 0
    # create one row per pair, where each row has 2 columns
    fig, ax = pl.subplots(len(pairs),2,figsize=(5.65,3), squeeze=False)
    for A, B in pairs:
        imA = ax[i][0].imshow(A, animated=True,vmin=0,cmap='Greys')
        imB = ax[i][1].imshow(B, animated=True,vmax=1,cmap='Greys')
        ax[i][0].axis('off')
        ax[i][1].axis('off')
        ax[i][0].set_title('A'+str(i))
        ax[i][1].set_title('B'+str(i))
        i += 1

    return fig, imA, imB


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

  
def simulate(genome, gridSize=200, N_simulation_steps=1000, delta_t=1.0):
    N = gridSize
    morphogens = [np.asarray(get_initial_A_and_B(N)) for i in range(len(genome))]
    
    for step in range(N_simulation_steps):
        # calculate parameters (important for dependancies)
        fValues = []
        kValues = []
        DAValues = []
        DBValues = []
        for chromosome in genome:
            fValues.append(decode(gene["f"], morphogens))
            kValues.append(decode(gene["k"], morphogens))
            DAValues.append(decode(gene["DA"], morphogens))
            DBValues.append(decode(gene["DB"], morphogens))
        
        
        for morphogen, f, k, DA, DB in zip(morphogens, fValues, kValues, DAValues, DBValues):
            morphogen = update(*morphogen, DA, DB, f, k, delta_t)
            
            
    draw( morphogens )

    # show the result
    pl.show()
    return morphogens
 
 
if __name__=="__main__":
    genome =                \
    [                       \
        defaultGeneSet(),   \
        defaultDependancyGeneSet() \
    ]
    simulate(genome)
