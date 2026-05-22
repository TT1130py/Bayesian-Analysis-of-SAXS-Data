#####----- Compare radius of gyration distributions between multiple sets of reweighed structures

#####----- IMPORTS
import numpy as np
import pandas as pd
import os
import argparse
import glob
from natsort import natsorted

#####----- PATHS, ARGUMENTS, REQUIREMENTS
experimental_data = "/users/t/j/tjaglal/experimental_data/SASDLU4.dat"
experimental_rg = 6.7 * 10 #convert to angstrom

## Edit here and/or add another subset for every weight file used
post_weights_1 = ""
grid_sum_1 = "" #this is grid sum
save_path_1 = ""
best_dro_1 = ""
best_r0_1 = ""

## Amount of weights used
wam = ""

#####----- FUNCTIONS

def concat_rg(best_dro, best_r0, grid_path, save_path):
    grid_df = pd.read_csv(grid_path, sep=',')
    gp = grid_df.index[(grid_df['d_rho'] == best_dro) & (grid_df['r0'] == best_r0)].tolist()[0]

    search_pattern = os.path.join(save_path, "mm*", f"GP{gp}", "Rg_env.dat")
    found_files =  glob.glob(search_pattern)

    sorted_files = natsorted(found_files)
    rg_list = []

    for file in sorted_files:
        data = pd.read_csv(file, sep='\s+', header=None, names=['Rg'])

        parts = file.split(os.sep)
        frac_folder = [p for p in parts if p.startswith('mm')][0]
        system, suffix = frac_folder.split('_')

        if suffix == "end":
            data['frame'] = range(len(data))
            data['is_end'] = True
        else:
            block_end = int(suffix)
            block_start = block_end - len(data) + 1
            data['frame'] = data.index + block_start
            data['is_end'] = False

        data['system'] = system
        rg_list.append(data)

    rg_df = pd.concat(rg_list, ignore_index=True)
    for system in rg_df['system'].unique():
        mask_sys = rg_df['system'] == system
        mask_end = mask_sys & (rg_df['is_end'])

        if mask_end.any():
            max_frame = rg_df.loc[mask_sys & (~rg_df['is_end']), 'frame'].max()
            n_end = mask_end.sum()

            rg_df.loc[mask_end, 'frame'] = range(int(max_frame) + 1, int(max_frame) + 1 + n_end)

    rg_df = rg_df.drop(columns=['is_end'])

    return rg_df