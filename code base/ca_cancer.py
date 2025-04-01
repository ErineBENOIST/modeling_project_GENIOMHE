from cellularautomata import CountType, GuiCA
from random import random

A0 = 0.1           # ATP threshold
PA = 0.001         # Mutation rate
TFIRE = 7          # Used for visual debugging, not fire here (I use fire forest as references)

def CA_Cancer(cell, neighbors: list):
    category, traits = cell
    g = 1.0  # glucose placeholder
    c = 1.0  # oxygen placeholder
    k = 5.0  # glycolytic boost
    n = 1/18

    glycolytic = 'G' in traits
    acid_resist = 'A' in traits
    hyperplastic = 'H' in traits

    # Metabolism
    delta_g = g * (k if glycolytic else 1)
    delta_a = c + n * (delta_g - c)
    delta_h = delta_g - c

    # Cell death by low ATP
    if delta_a < A0:
        return ('Empty', '')

    # Death by acidity (simplified threshold)
    if (not acid_resist and delta_h > 0.5) or (acid_resist and delta_h > 1.5):
        return ('Empty', '')

    # Division condition
    can_divide = any(nc[0] == 'Empty' for nc in neighbors)
    if can_divide and random() < ((delta_a - A0) / (1 - A0)):
        new_traits = traits
        for t in ['H', 'G', 'A']:
            if random() < PA:
                new_traits = new_traits.replace(t, '') if t in new_traits else new_traits + t
        return ('Cell', new_traits)

    # Default: stay the same
    return (category, traits)

# Colors per phenotype (not all combinations are used)
cellcolors = {
    ('Empty', None): 'white',
    ('Cell', ''): 'gray',
    ('Cell H', 'H'): 'pink',
    ('Cell G', 'G'): 'green',
    ('Cell A', 'A'): 'blue',
    ('Cell HG', 'HG'): 'lime',
    ('Cell GA', 'GA'): 'yellow',
    ('Cell HA', 'HA'): 'purple',
    ('Cell HGA', 'HGA'): 'orange',
}

GuiCA(CA_Cancer, cellcolors)
