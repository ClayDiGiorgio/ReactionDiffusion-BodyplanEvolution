###
# TODO List
###
# - add a field in the chromosome saying which index of 
#   morphogen it codes for.
#       * Right now, a haploid genome always codes for morphogens 0 to N-1
#         they can't skip morphogen 2, for example
# - add some way to linearize each chromosome
# - support reproduction between two diploid genomes, including crossover
#   and chromosome duplication / deletion


#from benmaier_reactionDiffusion.gray_scott_static import draw
from benmaier_reactionDiffusion.gray_scott_static import get_initial_A_and_B
from benmaier_reactionDiffusion.gray_scott_static import update
import matplotlib.pyplot as pl
import numpy as np
import random


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


# returns a diploid genome, eg
#   ( 
#    [                              # +-----< haploid genome 1 (genomeF)
#     {                           \ # | +--< Chromosome 1
#       "f": {"-":0.060},         \ # | |
#       "k": {"-":0.062},         \ # | |
#       "DA":{"-":0.16},          \ # | |
#       "DB":{"-":0.08}           \ # | |
#     }                             # | +-->
#    ],                             # +----->
#
#    [                              # +-----< haploid genome 2 (genomeM)
#     {                           \ # | +--< Chromosome 1
#       "f": {"-":0.060},         \ # | |
#       "k": {"-":0.062},         \ # | |
#       "DA":{"-":0.16},          \ # | |
#       "DB":{"-":0.08}           \ # | |
#     },                            # | +-->
#     {                           \ # | +--< Chromosome 2
#       "f": {"-":0.060},         \ # | |
#       "k": {"-":0.062},         \ # | |
#       "DA":{"-":0.16},          \ # | |
#       "DB":{"-":0.08}           \ # | |
#     }                             # | +-->
#    ]                              # +----->
#   )
def randomGenome(N):
    haploid1 = [defaultChromosome()]
    for i in range(N-1):
        haploid1.append(randomChromosome(N))
        
    haploid2 = [defaultChromosome()]
    for i in range(N-1):
        haploid2.append(randomChromosome(N))
    
    genome = (haploid1, haploid2)
    return genome


def simpleDefaultGenome(N):
    haploid1 = [defaultChromosome() for i in range(N)]
    haploid2 = [defaultChromosome() for i in range(N)]
    
    genome = (haploid1, haploid2)
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


# genomeF and genomeM are haploid genomes (they only contain at most one of each chromosome)  
def simulate(genomeF, genomeM, gridSize=200, N_simulation_steps=1000, delta_t=1.0):
    N = gridSize
    #morphogens = [np.asarray(get_initial_A_and_B(N)) for i in range(max(len(genomeF), len(genomeM)))]
    # morphogens is of the form [np.array([][], [][])]
    numGenomes = max(len(genomeF), len(genomeM))
    morphogens = [tuple(np.asarray(L) for L in get_initial_A_and_B(N)) for i in range(numGenomes)]
    # morphogens is of the form [(np.array(2d), np.array(2d)), ...]
    
    #print(morphogens)
    #print("=============")
    #print(morphogens[0])
    #print("=============")
    #print(morphogens[0][0].shape)
    #print(numGenomes)
    
    for step in range(N_simulation_steps):
        newMorphogens = []
        
        #print("\n===== Step " + str(step) + " ======")
        for haploid in (genomeF, genomeM):
            thisHaploidNewMorphogens = []
            # calculate parameters (important for calculating dependancies properly)
            fValues = []
            kValues = []
            DAValues = []
            DBValues = []
            for chromosomeF, chromosomeM in zip(genomeF, genomeM):
                fValues. append(decode(chromosomeF["f"],  morphogens) + decode(chromosomeM["f"],  morphogens))
                kValues. append(decode(chromosomeF["k"],  morphogens) + decode(chromosomeM["k"],  morphogens))
                DAValues.append(decode(chromosomeF["DA"], morphogens) + decode(chromosomeM["DA"], morphogens))
                DBValues.append(decode(chromosomeF["DB"], morphogens) + decode(chromosomeM["DB"], morphogens))
            
            #i = 0
            for morphogen, f, k, DA, DB in zip(morphogens, fValues, kValues, DAValues, DBValues):
                #print("\nChromosome/morphogen " + str(i))
                #i += 1
                
                A = morphogen[0]
                B = morphogen[1]
                f = np.clip(f, 0, 1)
                k = np.clip(k, 0, 1)
                DA = np.clip(DA, 0, 1)
                DB = np.clip(DB, 0, 1)
                
                #print("A = array" + str(A.shape))
                #print("B = array" + str(A.shape))
                #if type(f) == type(A):
                    #print("f = array" + str(f.shape))
                #if type(k) == type(A):
                    #print("k = array" + str(k.shape))
                #if type(DA) == type(A):
                    #print("DA = array" + str(DA.shape))
                #if type(DB) == type(A):
                    #print("DB = array" + str(DB.shape))
                
                A,B = update(A, B, DA, DB, f, k, delta_t)
                #print("\t" + str(A.shape))
                thisHaploidNewMorphogens.append((A, B))
                
            newMorphogens.append(thisHaploidNewMorphogens)
        
        # newMorphogens is of the form [[(A,B), ...], [(A, B), ...]]
        # where newMorphogens[0] is a list of all updated morphogen A, B pairs
        # updated by genomeF
        # and newMorphogens[1] is a list of all updated morphogens by genomeM
        
        # morphogens = [(0.5*newFPair[0]+0.5*newMPair[0], 0.5*newFPair[1] + 0.5*newMPair[1]) for newFPair, newMPair in zip(newMorphogens[0], newMorphogens[1])]
        morphogens = newMorphogens[0]
        print(morphogens)
            
    draw( morphogens )

    # show the result
    pl.show()
    return morphogens
 

if __name__=="__main__":
    #genome =                \
    #[                       \
        #defaultGeneSet(),   \
        #defaultDependancyGeneSet() \
    #]
    #simulate(genome)
    
    simulate(*simpleDefaultGenome(4))
    
    #import pprint
    #pp = pprint.PrettyPrinter(indent=4)
    
    #genome = randomGenome(4)
    #pp.pprint(genome)
    #simulate(*genome)
