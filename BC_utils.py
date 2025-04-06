
from random import random, choice, seed
seed(10)


# constants
dg = 1.3e2
dc = 5
k = 10

a0 = 0.1    # ATP production threshold, if lower than alpha0, cell dies
hN = 9.3e2   # threshold of local H+ level for normal cells, die if greater
hT = 8.6e3  # threshold of local H+ level for acid-resistant cells, die if greater
pa = 1e-3    # probability of randomly acquiring one of the three phenotypes (H,G,A)

moores = [
        (-1,-1), (-1,0), (-1,1),
        (0,-1),          (0,1),
        (1,-1),  (1,0), (1,1)
    ]


def UpdateMetabolites(phenotype, neighbors):
    
    # updating glucose level
    if phenotype == "empty":
        deltaG = 0 # all delta values = 0 in vacant cells
    elif "G" in phenotype:
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
        _, env = neighbor
        sum_gluc += env[0]
        sum_oxy += env[1]
        sum_acid += env[2]

    gluc_level = sum_gluc / (4 + deltaG)
    oxy_level = sum_oxy / (4 + deltaO)

    
    # updating H+ level
    if phenotype == "empty":
        deltaH = 0
    elif "G" in phenotype:
        deltaH = k * gluc_level - oxy_level
    else:
        if gluc_level > oxy_level:
            deltaH = gluc_level - oxy_level
        else:
            deltaH = 0  
            
    h_level = (sum_acid + deltaH) / 4

    return gluc_level, oxy_level, h_level



def select_daughter_neighbor(oxygen_in_neighbors):
    """
    Selecting the neighbor with the highest oxygen level for the daughter cell
    If there are more than one cell (n) with the highest oxygen level, choose one among them randomly
    """    
    max_oxygen = max(oxygen_in_neighbors.values()) 
    max_oxygen_indices = [k for k, v in oxygen_in_neighbors.items() if v == max_oxygen]  # get all keys with the max oxygen level
    
    return moores[choice(max_oxygen_indices)] # random.choice


def acquire_phenotypes0(parent_phenotype, p_a=pa):
    """
    Applies phenotype changes to the daughter cells. Each daughter cell has a probability of p_a 
    to acquire one of new traits (A, G, H) or losing (A,G,H) because changes are reversible
    Changes can be:
      1) Gaining one new trait A,G,H if parent is normal
      2) Losing one existing trait (e.g., AH->H, AGH->AG, G->normal, etc.)
      3) Switching one existing trait for another (e.g. A->G, G->H, etc.)
    """
    assert parent_phenotype != "empty"
    
    traits = set(parent_phenotype.replace("normal", ""))
    all_traits = {"A", "G", "H"}

    if random() < p_a:
        
        # parent is normal (no trait = "")
        if not traits:
            traits.add(choice(list(all_traits)))
    
        # parent is AGH
        elif len(traits) == 3:
            traits.remove(choice(list(traits)))
    
        # parent is A, G, H, AG, AH, GH
        else:
            action = choice(["add", "remove", "swap"])
            # swap: A->G, G->H, etc.
            # add: A->AG, G->AG, GH->AGH etc.
            # remove: A->normal, GH->G
    
            if action == "remove":
                traits.remove(choice(list(traits)))
            elif action == "swap":
                old_trait = choice(list(traits))
                new_trait = choice(list(all_traits - traits))
                traits.remove(old_trait)
                traits.add(new_trait)
            else:
                traits.add(choice(list(all_traits - traits)))
            
    return ''.join(sorted(traits)) if traits else "normal"

    
def acquire_phenotypes2(parent_phenotype, p_a=pa):
    """
    Strictly follows paper's description:
    - Each daughter cell INDEPENDENTLY has p_a chance to toggle ONE trait
    - Toggle = add if absent, remove if present
    - All traits (A/G/H) have equal probability of being selected
    """
    traits = set(parent_phenotype.replace("normal", ""))
    all_traits = {'A', 'G', 'H'}
    
    if random() < p_a:
        # Randomly select any one trait (A/G/H)
        selected_trait = choice(list(all_traits))
        
        # Toggle: Add if absent, Remove if present
        if selected_trait in traits:
            traits.remove(selected_trait)
        else:
            traits.add(selected_trait)
    
    return ''.join(sorted(traits)) if traits else "normal"


def acquire_phenotypes(parent_phenotype, p_a=pa):
    """
    Applies phenotype changes to the daughter cells. Each daughter cell has a probability of p_a 
    to acquire one of new traits (A, G, H) or losing (A,G,H) because changes are reversible.
    Changes can be:
      1) Gaining one new trait A,G,H if parent is normal
      2) Losing one existing trait (e.g., AH->H, AGH->AG, G->normal, etc.)
      3) Switching one existing trait for another (e.g. A->G, G->H, etc.)
    """
    traits = set(parent_phenotype.replace("normal", ""))
    all_traits = {'A', 'G', 'H'}
    
    if random() < p_a:
        if not traits:
            # If the cell is "normal", it can only gain a trait
            new_trait = choice(list(all_traits))
            traits.add(new_trait)
        else:
            # If the cell has traits, decide whether to gain, lose, or switch
            action = choice(['gain', 'lose', 'switch'])
            
            if action == 'gain':
                # Gain a new trait not currently present
                possible_gains = all_traits - traits
                if possible_gains:
                    new_trait = choice(list(possible_gains))
                    traits.add(new_trait)
            
            elif action == 'lose':
                # Lose an existing trait
                if traits:
                    trait_to_remove = choice(list(traits))
                    traits.remove(trait_to_remove)
            
            elif action == 'switch':
                # Switch one existing trait for another
                if traits:
                    trait_to_switch = choice(list(traits))
                    traits.remove(trait_to_switch)
                    possible_switches = all_traits - {trait_to_switch} - traits
                    if possible_switches:
                        new_trait = choice(list(possible_switches))
                        traits.add(new_trait)
    
    return ''.join(sorted(traits)) if traits else "normal"


def get_targeting_neighbor(neighbors):
    """
    Check if current cell (with target=(None, None)) is targeted by any neighbors.
    If multiple neighbors target this cell, randomly select one.
    Returns:
        The neighbor cell object that targets this cell, or None if none do
    """
    targeting_neighbors = []
    
    # Moore neighborhood relative positions in the specified order
    
    for idx, neighbor in enumerate(neighbors):
        # Extract target coordinates (the [e,f] list)
        target_coords = neighbor[1][3][0]  # Gets the [e,f] list
        
        # Skip if neighbor has no target (shouldn't happen with your structure)
        if target_coords is None:
            continue
        else:    
            # Get this neighbor's relative position to current cell
            if isinstance(target_coords, int):
                #print(target_coords)
                break
            required_target = [a + b for a, b in zip(moores[idx], target_coords)]
            
            if required_target == [0,0]:
                targeting_neighbors.append(neighbor)
    
    if not targeting_neighbors:
        return None
    
    return choice(targeting_neighbors)


