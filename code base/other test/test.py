# This code is without the prof interface
from random import random, choice
import numpy as np
from scipy.ndimage import convolve
import matplotlib.pyplot as plt

# config

GRID_WIDTH = 50  # N
CELL_DIAMETER = 25e-6  # 25 µm
TIME_STEPS = 200
A0 = 0.1  # ATP threshold for survival
PA = 1e-3  # mutation probability

# Diffusion coefficients
DG = 5e-6  # Glucose (cm²/s)
DC = 1.46e-5  # Oxygen (cm²/s)
DH = 1.08e-5  # Acid proxy (cm²/s)

# Consumption rates
K_NG = 5e-5
K_TG = 5e-4
K_C = 9.41e-2

# Non-dimensional constants
N = 1 / 18
K = K_TG / K_NG

# cell

class Cell:
    def __init__(self, attached=True):
        self.attached = attached
        self.proliferative = False
        self.glycolytic = False
        self.acid_resistant = False

    def mutate(self):
        traits = ['proliferative', 'glycolytic', 'acid_resistant']
        if random() < PA:
            trait = choice(traits)
            setattr(self, trait, not getattr(self, trait))

    def can_survive(self, atp, h_conc):
        if atp < A0:
            return False
        if not self.acid_resistant and h_conc > 9.3e2:
            return random() > 9.3e2 / h_conc
        elif self.acid_resistant and h_conc > 8.6e3:
            return random() > 8.6e3 / h_conc
        if not self.attached and not self.proliferative:
            return False
        return True

    def divide_prob(self, atp):
        if atp < A0:
            return 0
        return min(1, (atp - A0) / (1 - A0))
    

# grid

class Grid:
    def __init__(self):
        self.grid = [[Cell() if i == 0 else None for _ in range(GRID_WIDTH)] for i in range(GRID_WIDTH)]
        self.oxygen = np.ones((GRID_WIDTH, GRID_WIDTH))
        self.glucose = np.ones((GRID_WIDTH, GRID_WIDTH))
        self.acid = np.zeros((GRID_WIDTH, GRID_WIDTH))

    def neighbors(self, x, y):
        offsets = [(-1,0), (1,0), (0,-1), (0,1)]
        return [(x+dx, y+dy) for dx, dy in offsets if 0 <= x+dx < GRID_WIDTH and 0 <= y+dy < GRID_WIDTH]

    def step(self, update_callback):
        new_grid = [[None for _ in range(GRID_WIDTH)] for _ in range(GRID_WIDTH)]
        for i in range(GRID_WIDTH):
            for j in range(GRID_WIDTH):
                cell = self.grid[i][j]
                if cell:
                    c = self.oxygen[i][j]
                    g = self.glucose[i][j]
                    delta_g = g * (5 if cell.glycolytic else 1)
                    atp = c + (1/18) * (delta_g - c)
                    acid = delta_g - c
                    h_conc = acid * 1000  # proxy
                    if not cell.can_survive(atp, h_conc):
                        continue
                    if random := cell.divide_prob(atp) > np.random.rand():
                        for ni, nj in self.neighbors(i, j):
                            if not self.grid[ni][nj] and not new_grid[ni][nj]:
                                new_cell = Cell(attached=False)
                                new_cell.proliferative = cell.proliferative
                                new_cell.glycolytic = cell.glycolytic
                                new_cell.acid_resistant = cell.acid_resistant
                                new_cell.mutate()
                                new_grid[ni][nj] = new_cell
                                break
                    new_grid[i][j] = cell
        self.grid = new_grid
        update_callback(self)

# diffusion.py

def update_diffusion(grid):
    G = np.copy(grid.glucose)
    C = np.copy(grid.oxygen)
    H = np.copy(grid.acid)

    kernel = np.array([[0,1,0],[1,-4,1],[0,1,0]])

    dG = DG / (K_NG * CELL_DIAMETER**2)
    dC = DC / (K_NG * CELL_DIAMETER**2)
    dH = DH / (K_NG * CELL_DIAMETER**2)

    G += dG * convolve(G, kernel, mode='nearest')
    C += dC * convolve(C, kernel, mode='nearest')
    H += dH * convolve(H, kernel, mode='nearest')

    grid.glucose = np.clip(G, 0, 1)
    grid.oxygen = np.clip(C, 0, 1)
    grid.acid = np.clip(H, 0, 1)

# main.py

def get_grid_state(grid):
    state = np.zeros((GRID_WIDTH, GRID_WIDTH))
    for i in range(GRID_WIDTH):
        for j in range(GRID_WIDTH):
            cell = grid.grid[i][j]
            if cell:
                if cell.glycolytic and cell.acid_resistant:
                    state[i, j] = 4  # yellow
                elif cell.glycolytic:
                    state[i, j] = 3  # green
                elif cell.acid_resistant:
                    state[i, j] = 2  # blue
                elif cell.proliferative:
                    state[i, j] = 1  # pink
                else:
                    state[i, j] = 5  # normal (grey)
    return state

def simulate():
    g = Grid()
    fig, ax = plt.subplots()
    cmap = plt.get_cmap('tab10', 6)
    im = ax.imshow(get_grid_state(g), cmap=cmap, vmin=0, vmax=5)
    ax.set_title("Cellular Automaton Simulation")
    cbar = plt.colorbar(im, ax=ax, ticks=range(6))
    cbar.ax.set_yticklabels([
        'Empty', 'Proliferative', 'Acid-Resistant',
        'Glycolytic', 'Gly+AcidRes', 'Normal'
    ])

    for step in range(TIME_STEPS):
        g.step(update_diffusion)
        im.set_data(get_grid_state(g))
        ax.set_xlabel(f"Time step: {step+1}")
        plt.pause(0.1)

    plt.show()

if __name__ == "__main__":
    simulate()



