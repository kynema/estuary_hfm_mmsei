import numpy as np

n_cell = np.array([320, 256, 32]) # Base case
#n_cell = 2*n_cell
n_cell = 3*n_cell
max_level = 2           
blocking_factor_z = 4

# Domain 
prob_lo =  np.array([0, 0,  -170])
prob_hi =  np.array([6400, 5000, 42])

# Check that n_cell[2] is divisible by blocking_factor_z
if n_cell[2] % blocking_factor_z != 0:
    raise ValueError("n_cell[2] must be divisible by blocking_factor_z")

# Base resolution
lengths = prob_hi - prob_lo
base_resolution = lengths / n_cell

# Finest resolution at max_level
finest_resolution = base_resolution / (2 ** max_level)

# Print the results
print("n_cell:", n_cell)
print("Number of cells (level 0):", n_cell[0] * n_cell[1] * n_cell[2])
print("Base resolution (level 0):", base_resolution)
print("Finest resolution (level {}):".format(max_level), finest_resolution)