
from cellularautomata import CountType, GuiCA
from cellularautomata import GenerateCA, SimulateCA, ShowSimulation
from random import random, choice


### ---------- Rules for the metabolites ----------
# constants
dg = 1.3e2
dc = 5
k = 10

a0 = 0.1        # ATP production threshold, if lower than alpha0, cell dies
hN = 9.3e2          # threshold of local H+ level for normal cells, die if greater
hT = 8.6e3          # threshold of local H+ level for acid-resistant cells, die if greater
pa = 1e-3           # probability of randomly acquiring one of the three phenotypes (H,G,A)


def Metabolites(phenotype, neighbors):
    
    # updating glucose level
    if phenotype == "empty":
        deltaG = 0 # all delta values = 0 in vacant cells
    elif phenotype == "G":
        deltaG = k/(dg**2)
    else: # normal / non-glycolytic cells
        deltaG = 1/(dg**2)

    
    # updating oxygen level
    if phenotype == "empty":
        deltaO = 0
    else:
        deltaO = 1/(dc**2)


    sum_gluc = 0
    sum_oxy = 0
    sum_acid = 0
    for i in [1,3,4,6]:
        neighbor = neighbors[i] # the default neighbors is of Moore type -> get VonNeumman type using indices
        _, env in neighbor
        sum_gluc += env[0]
        sum_oxy += env[1]
        sum_acid += env[2]

    gluc_level = sum_gluc / (4+deltaG)
    oxy_level = sum_oxy / (4+deltaO)

    
    # updating H+ level
    if phenotype == "empty":
        deltaH = 0
    elif phenotype == "G":
        deltaH = gluc_level - oxy_level
    else:
        if gluc_level > oxy_level:
            deltaH = gluc_level - oxy_level
        else:
            deltaH = 0  
            
    h_level = sum_acid / (4+deltaH)

    return gluc_level, oxy_level, h_level


def select_daughter_neighbor(oxygen_in_neighbors):
    """
    Selecting the neighbor with the highest oxygen level for the daughter cell
    If there are more than one cell (n) with the highest oxygen level, choose one among them randomly
    """    
    max_oxygen = max(oxygen_in_neighbors.values()) 
    max_oxygen_indices = [k for k, v in oxygen_in_neighbors.items() if v == max_oxygen]  # get all keys with the max oxygen level
    
    return choice(max_oxygen_indices) # random.choice


def BC(cell, neighbors):
    phenotype, env = cell
    # cell = (phenotype, (gluc_level, oxy_level, h_level, daughter_index))
    # gluc_level, oxy_level, h_level = env # three elements of the environment are glucose, oxygen, and H+

    
    # ------------------ updating local levels of glucose, oxygen, and H+
    gluc_level, oxy_level, h_level = Metabolites(phenotype, neighbors)


    
    # ------------------ updating cell phenotype

    # if cell is empty and is the daughter index of its neighbor, then turn into that cell
    

    
    # ---------- cell death:
    if oxy_level < 0.05 and phenotype != "empty":
        return ("empty", (gluc_level, oxy_level, h_level))
    else:
        if phenotype != "A" and h_level < hN:
            p_death = h_level / hN
        elif phenotype == "A" and h_level < hT:
            p_death = h_level / hT
        else:
            p_death = 1

    die = random() < p_death
    if die:
        return ("empty", (gluc_level, oxy_level, h_level))
    # elif not die and :
    #     return ("A", (gluc_level, oxy_level, h_level))
        

    # ---------- cell division:
    # ----- ATP production:
    if phenotype == "G":
        phiG = k * gluc_level
    else:
        phiG = gluc_level
    phiA = oxy_level + 1/18 * (phiG - oxy_level)

    # ----- division probability:
    if phiA < a0:
        phenotype = "empty"
    elif phiA < 1 and phiA > a0:
        p_division = (phiA - a0) / (1-a0)
    elif phiA >=1:
        p_division = 1
    
    if random() < p_division:
        empty = CountType(neighbors, 'empty')

        if empty == 1:
            """the empty cell will be replaced by a new normal cell with p_div probability"""
            print("")
        elif empty > 1:
            oxygen_in_neighbors = {}
            for i in range(len(neighbors)):
                _, env = neighbors[i]
                if _ == "empty":
                    oxygen_in_neighbors[i] = env[1]

            # selecting the neighbor with highest oxygen level for the daughter cell
            daughter_index = select_daughter_neighbor(oxygen_in_neighbors)
            """one cell with the highest oxygen level will be replaced by a new normal cell with p_div probability"""
            print("")
        else:
            """the cell will fall into quiescent state"""
            return (phenotype, (gluc_level, oxy_level, h_level))
    else:
        return (phenotype, (gluc_level, oxy_level, h_level))


    ####### 1. complete division probability and occupation of new cells in the neighboring area
    ####### 2. switching phenotype according to environmental factors
            

    return (phenotype, (gluc_level, oxy_level, h_level))