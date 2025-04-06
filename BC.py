from BC_utils import (UpdateMetabolites, select_daughter_neighbor, 
                      acquire_phenotypes, get_targeting_neighbor)
from cellularautomata_BC import GuiCA
from random import random, seed

seed(10)

# constants
dg = 1.3e2
dc = 5
k = 10

a0 = 0.1        # ATP production threshold, if lower than alpha0, cell dies
hN = 9.3e2      # threshold of local H+ level for normal cells, die if greater
hT = 8.6e3      # threshold of local H+ level for acid-resistant cells, die if greater
pa = 1e-3      # probability of randomly acquiring one of the three phenotypes (H,G,A)

moores = [
        (-1,-1), (-1,0), (-1,1),
        (0,-1),          (0,1),
        (1,-1),  (1,0), (1,1)
    ]


def BC(cell, neighbors):
    phenotype, env = cell
    
    # ------------------ 1. UPDATING LEVELS OF GLUCOSE, O2, AND H+ ------------------
    gluc_level, oxy_level, h_level = UpdateMetabolites(phenotype, neighbors)
    
    # ------------------ updating empty element
    if phenotype == "empty":
        chosen_by = get_targeting_neighbor(neighbors)
        if chosen_by is None: 
            return (phenotype, (gluc_level, oxy_level, h_level, (None, None)))
        else:
            new_phenotype = chosen_by[1][3][1]
            return (new_phenotype, (gluc_level, oxy_level, h_level, (None, None)))
            
    # ------------------ 2. CELL DEATH ------------------
    if "A" in phenotype:
        h_threshold = hT
    else:
        h_threshold = hN
    
    if h_level < h_threshold:
        p_death = h_level / h_threshold
    else:
        p_death = 1
    
    if random() < p_death:
        return ("empty", (gluc_level, oxy_level, h_level, (None, None)))

    # ------------------ 4. CELL DIVISION ------------------
    if "G" in phenotype:
        phiG = k * gluc_level
    else:
        phiG = gluc_level
    phiA = oxy_level + (phiG - oxy_level) / 18

    # cell will die if produce ATP (phiA) < a0
    if phiA < a0:
        return ("empty", (gluc_level, oxy_level, h_level, (None, None)))
    elif phiA < 1 and phiA > a0:
        p_division = (phiA - a0) / (1-a0)
    elif phiA >=1:
        p_division = 1

    if not random() < p_division: # no division, stay quiescent (same phenotype)
        return (phenotype, (gluc_level, oxy_level, h_level, (None, None)))
        
    else:
        empty_neighbors_o2 = {}
        for i in range(len(neighbors)):
            _, env = neighbors[i]
            if _ == "empty":
                empty_neighbors_o2[i] = env[1] # store oxygen level of empty neighbors

        if len(empty_neighbors_o2) == 0:
            return (phenotype, (gluc_level, oxy_level, h_level, (None, None)))
        elif len(empty_neighbors_o2) == 1:
            daughter_index = list(empty_neighbors_o2.keys())[0]
        else: # if more than 2 empty neighbor exists, choose one with highest O2
            daughter_index = select_daughter_neighbor(empty_neighbors_o2)

        if daughter_index is not None: 
            # if a location is found for daughter cells, choose phenotype
            daughter1_phenotype = acquire_phenotypes(phenotype)
            daughter2_phenotype = acquire_phenotypes(phenotype)
            return (daughter1_phenotype, (gluc_level, oxy_level, h_level, (daughter_index, daughter2_phenotype)))
        else:
            return (phenotype, (gluc_level, oxy_level, h_level, (None, None)))

# ===================== Main program =====================
N = 50
M = 50
g_base = 1.0 
c_base = 1.0
h_base = 0.0
# cellcolors = {('empty', (0.0, 0.0, 0.0, (None, None))): 'white', 
#               ('normal', (g_base, g_base, h_base, (None, None))): 'grey', 
#               ('H', (None, None, None, (None, None))): '#f38e8d',    # pink
#               ('G', (None, None, None, (None, None))): '#09a44c',    # green
#               ('GH', (None, None, None, (None, None))): '#09a44c',   # green
#               ('A', (None, None, None, (None, None))): '#56529e',    # blue
#               ('AH', (None, None, None, (None, None))): '#56529e',   # blue
#               ('AG', (None, None, None, (None, None))): 'yellow', 
#               ('AGH', (None, None, None, (None, None))): 'black'} 

# cellcolors = {('empty', (0.0, 0.0, 0.0, (None, None))): 'white', 
#               ('normal', (g_base, g_base, h_base, (None, None))): 'grey', 
#               ('H', (None, None, None, (None, None))): '#f38e8d',    # pink
#               ('G', (None, None, None, (None, None))): '#3ab670',    # green
#               ('GH', (None, None, None, (None, None))): '#067335',   # dark green
#               ('A', (None, None, None, (None, None))): '#7875b1',    # blue
#               ('AH', (None, None, None, (None, None))): '#45427e',   # dark blue
#               ('AG', (None, None, None, (None, None))): 'yellow', 
#               ('AGH', (None, None, None, (None, None))): 'black'}

# cellcolors = {('empty', (0.0, 0.0, 0.0, (None, None))): 'white', 
#               ('normal', (g_base, g_base, h_base, (None, None))): 'grey', 
#               ('H', (None, None, None, (None, None))): '#f38e8d',    # pink
#               ('G', (None, None, None, (None, None))): '#3ab670',    # green
#               ('GH', (None, None, None, (None, None))): 'red',   # dark green
#               ('A', (None, None, None, (None, None))): 'blue',    # blue
#               ('AH', (None, None, None, (None, None))): '#45427e',   # dark blue
#               ('AG', (None, None, None, (None, None))): 'yellow', 
#               ('AGH', (None, None, None, (None, None))): 'black'}

cellcolors = {('empty', (0.0, 0.0, 0.0, (None, None))): 'white', 
              ('normal', (g_base, g_base, h_base, (None, None))): 'grey', 
              ('H', (None, None, None, (None, None))): '#f38e8d',    # pink
              ('G', (None, None, None, (None, None))): 'black',   
              ('GH', (None, None, None, (None, None))): '#067335',   # dark green
              ('A', (None, None, None, (None, None))): '#45427e',    # blue
              ('AH', (None, None, None, (None, None))): '#45427e',   # dark blue
              ('AG', (None, None, None, (None, None))): 'black', 
              ('AGH', (None, None, None, (None, None))): '#eecb4a'} # yellow

GuiCA(BC, cellcolors, gridsize=100, duration=800)
