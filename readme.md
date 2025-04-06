# Cellular automata for the somatic evolution of breast cancer cells

**M1 GENIOMHE - Univ Evry Paris-Saclay (2024-2025)**

- Minh Ngoc VU

- Erine Benoist

This github repo is a supplementary material for the our final report for our course ***Modeling and simulation of biological system***. We reproduced the cellular automaton described in the paper "Cellular adaptations to hypoxia and acidosis during somatic evolution of breast cancer" (Gatenby, Robert A., et al., British journal of cancer 97.5 (2007): 646-653).

This repo contains scripts for simulations, data, and figures that we included in our report. The framework used for the simulation was adapted from our Professor Delaplace's [Cellular-Automata-Python library](https://github.com/Franck-Delaplace/Cellular-Automata-Python). 


### Requirements
The simulation was conducted on `python=3.10.8` and `numpy=1.23.5`. It will fail if newer python or numpy versions are used.

### Scripts
Three main scripts were used for the simulations:
- `BC.py`: our main script. The function `BC` describe the rules for cell dynamics that should be applied for each cell in the automaton
- `BC_utils.py`: utility script, containing helper functions for: 1) updating the metabolite level in each cell (glucose, oxygen, and H+), 2) phenotype acquisition for daughter cells during division, and 3) selecting neighbor destination for daughter cell placement
- `cellularautomata_BC.py`: adapted from the original library's `cellularautomata.py`. The GenerateCA_BC and SimulationCA_BC were created to handle row-specific rules for the CA, for example, dealing with the basement membrane (bottom layer of the grid). Additional code was made to save the cell count data from the simulation to .csv files. Other modifications concern plots and fonts. 

### Rules:
Detailed rules of the cellular automaton and the source paper can be found in [ref](ref/). 

### Simulation results
We successfully found 2 evolution pathways described in the model. In addition, we also identified a novel pathway (called X) that shows the co-evolution of three invasive phenotypes in the cell system. Simulation xamples can be found in [gif](gif/).

**Pathway 1:**

[pathway1](gif/pathway1.gif)

**Pathway 2:**

[pathway2](gif/pathway2.gif)

**Pathway X (new, not yet validated):**

[pathwayX](gif/pathwayX2_a0.05.gif)