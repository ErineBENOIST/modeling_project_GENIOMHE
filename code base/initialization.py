from cellularautomata import GenerateCA, GuiCA
from ca_cancer import CA_Cancer
import numpy as np

# Cell colors and categories
cellcolors = {
    ('Empty', ''): 'white',
    ('Cell', ''): 'gray',         # Normal cells will be 'gray'
}

# Create the grid with a size of 50x50 (adjust as needed)
n = 50  # Grid size (NxN)

# Create the grid and set the first row to be normal cells
initial_grid = np.array([[('Cell', '') if j == 0 else ('Empty', '') for i in range(n)] for j in range(n)])

# Pass the initial grid to GuiCA and assign the animation to a variable
anim = GuiCA(CA_Cancer, cellcolors, initial_grid=initial_grid)
