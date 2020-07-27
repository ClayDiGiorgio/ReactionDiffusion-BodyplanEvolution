
# should be a list of elements, where each element is one of the following
# }
# {
# ~
# {A|B}{0-9}+
# a floating point number
#

def chromosomeFromLinearization(linearChromo):
    chromosome = {}
    genes = ["f", "k", "DA", "DB"]
    
    bracketsOpen = False
    currentGene = {}
    i = -1
    for _ in range(len(linearChromo)):
        i += 1
        
        if not bracketsOpen:
            if linearChromo[i] == '{':
                bracketsOpen = True
            continue
        
        if linearChromo[i] == '}':
            bracketsOpen = False
            chromosome[genes.pop(0)] = currentGene
            currentGene = {}
            continue
        
        if linearChromo[i] == '~':
            if i == len(linearChromo)-1:
                continue
            i += 1
            currentGene['~'] = linearChromo[i]
            continue

        # the only case left is A0
        if type(linearChromo[i]) == type('string'): 
            if i == len(linearChromo)-1:
                continue
            currentGene[linearChromo[i]] = linearChromo[i+1]
            i += 1
