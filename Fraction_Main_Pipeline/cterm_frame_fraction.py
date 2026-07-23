############---------- REQUIREMENTS
import os
import subprocess
import sys
import pandas as pd
import numpy as np
from datetime import date
import time
import glob
import re

#Date Time
epoch = time.time()
cd = time.strftime("%Y%m%d-%H%M%S", time.localtime(epoch))
today = date.today()

############---------- PATHS
path_frames = "/users/t/j/tjaglal/cterm_structures"
folder_pattern = os.path.join(path_frames, "MD*", "MD*_*")
structure_folders = glob.glob(folder_pattern)

#Optional filter
#ignore_system = ["mm128", "mm256"]
#structure_foldders = [folder for folder in structure_folders if not any (sys in folder for sys in ignore_system)]

path_exp_file = "/users/t/j/tjaglal/experimental_data"
path_main = "/users/t/j/tjaglal/Projects/iBME"
path_shell = "/users/t/j/tjaglal/shell"


############---------- FUNCTIONS

def create_grid(dro_min, dro_max, dro_step, r0_min, r0_max, r0_step):
    dro_grid = np.arange(dro_min, dro_max + dro_step, dro_step)
    r0_grid = np.arange(r0_min, r0_max + r0_step, r0_step)

    grid_data = []
    rm = 1.60

    index = 0
    for d in dro_grid:
        dro_val = round(d, 2)
        for r in r0_grid:
            r0_ratio = r / rm
            r0_val = round(r0_ratio, 4)

            grid_data.append([index, dro_val, r0_val])
            index += 1

    grid_name = np.array(grid_data)

    return grid_name, index, dro_min, dro_max, r0_min, r0_max

#############----------- MAIN

theta_val = 4000
gn, idx, d_min, d_max, r_min, r_max = create_grid(-26.72, 70.14, 3.34, 1.36, 1.76, .04)
gdf = pd.DataFrame(gn, columns=['index', 'dro', 'r0'])

#Make output file for this specific grid
path_results = "/gpfs1/home/t/j/tjaglal/Projects/iBME/cterm_grid_sims"
run_fol = "saxs_dro_{}_to_{}_r0_{}_to_{}".format(d_min, d_max, r_min, r_max)
sub_path = os.path.join(path_results, run_fol)
os.mkdir(sub_path)

np.savetxt(os.path.join(sub_path, "grid_full.txt"), gn, fmt=['%d', '%.2f', '%.2f'])
os.chdir(path_shell)

#Run Pepsi for fraction
fraction_job_ids = []

for struct_path in structure_folders:
    if os.path.isdir(struct_path):

        folder_name = os.path.basename(struct_path)
        spec_save = os.path.join(sub_path, folder_name)
        os.makedirs(spec_save, exist_ok=True)

        np.savetxt(
            os.path.join(spec_save, "grid_full.txt"),
            gn,
            fmt=['%d', '%.2f', '%.2f'])

        cmd = ["sbatch", "--parsable", "cterm_iBME_scan.sh", str(theta_val), spec_save, struct_path]
        run = subprocess.run(cmd, capture_output=True, text=True)

        if run.returncode == 0:
            job_id = run.stdout.strip()
            fraction_job_ids.append(job_id)
        else:
            print("Failed fulder submit")
    else:
        print("Pepsi Job failed")

