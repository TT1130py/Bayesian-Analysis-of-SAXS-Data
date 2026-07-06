import os
import glob
import pandas as pd
import numpy as np
import re
from natsort import natsorted
from datetime import date
import argparse

# --- Setup Arguments ---
parser = argparse.ArgumentParser(description="Recover and save structure weights")
parser.add_argument("save_path", type=str, help="Path to your save directory")
args = parser.parse_args()

# --- Reconstruct Paths ---
struc_path = "/users/t/j/tjaglal/structures"  
ibme_out_dir = "/gpfs1/home/t/j/tjaglal/Projects/iBME/iBME_res/iBME_dro_-26.72_to_13.36_r0r_1.36_to_1.68_t100"
grid_file_path = os.path.join(args.save_path, "grid_full.txt")
grid_sum_path = os.path.join(ibme_out_dir, "GRID_sum.txt")
today = date.today()

# --- Reload Data ---
print("Loading previously computed grid data...")
GRID_DF = pd.read_csv(grid_file_path, sep=r'\s+', header=None, names=['index', 'dro', 'r0'])
grid = np.genfromtxt(grid_sum_path, skip_header=1, delimiter=',', filling_values=np.nan) 

# --- Recalculate Best dro and r0 ---
chi2 = np.clip(grid[:,4], 1e-12, None)
phi  = np.clip(grid[:,5], 1e-12, None)
gamma = np.log(chi2 / phi)

# Find the exact row index of the minimum gamma
best_idx = np.nanargmin(gamma)
best_dro = grid[best_idx, 1]
best_r0 = grid[best_idx, 2]

print(f"Recovered Best Parameters -> δρ={best_dro:.2f}, r0={best_r0:.3f}")

# --- Execute Original SAVE Logic ---
weight_idx = GRID_DF.index[(GRID_DF['dro'] == best_dro) & (GRID_DF['r0'] == best_r0)].tolist()[0]
best_gp_dir = os.path.join(ibme_out_dir, f"GP{weight_idx}")

# Dynamically find the last .weights.dat file
weight_files = glob.glob(os.path.join(best_gp_dir, "*.weights.dat"))
if not weight_files:
    raise FileNotFoundError(f"No .weights.dat files found in {best_gp_dir}")
    
weight_files_sorted = sorted(weight_files, key=lambda x: int(re.search(r"_(\d+)\.weights\.dat", os.path.basename(x)).group(1)))
best_weight_file = weight_files_sorted[-1]

# Get a sorted list of ALL structure names to map the weights back to the PDBs
all_structures = glob.glob(os.path.join(struc_path, "*.pdb"))
contents = pd.DataFrame(natsorted([os.path.basename(x) for x in all_structures]))

# Map and save
opt_weight = pd.read_csv(best_weight_file, sep=r'\s+', header=None)
opt_weight['PDB_Name'] = opt_weight.iloc[:, 0].map(contents.iloc[:, 0])
opt_sorted = opt_weight.sort_values(by=1, ascending=False)

weights_out = os.path.join(ibme_out_dir, f'structure_weights_sorted_{today}.txt')
opt_sorted.to_csv(weights_out, index=None, sep='\t')

print(f"Success! Top structure weights saved to: {weights_out}")
