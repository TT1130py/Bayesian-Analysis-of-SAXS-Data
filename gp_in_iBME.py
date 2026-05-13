#gp in iBME
#running do_gp.sh from inside a python script rather than iBME in shell script
#RUN FROM VENV

import subprocess
import os
from gp_files import iBME_script
import pandas as pd
import glob 
import numpy as np
import re
import matplotlib.pyplot as plt
import time
from datetime import date
from natsort import natsorted
import shutil

##Create main and sub directory

#Time and date of run
epoch = time.time()
cd = time.strftime("%a, %d, %b, %Y, %H:%M:%S", time.localtime(epoch))
today = date.today()

main_path = '/home/malab/iBME'

sub_dir_num = 0
run_fol = "run_{}_{}"
sub_path = os.path.join(main_path, run_fol)
while os.path.isdir(sub_path.format(today, sub_dir_num)):
    sub_dir_num += 1
sub_path = sub_path.format(today, sub_dir_num)
os.mkdir(sub_path)

##Set paths and parameters

#Working directory
print("Current working directory: {0}".format(os.getcwd()))
os.chdir('/home/malab/iBME')
print("Current working directory: {0}".format(os.getcwd()))
os.chdir(sub_path)
print("Current working directory: {0}".format(os.getcwd()))

#Parameters for do_gp
path_structures = "/home/malab/Desktop/PDB_ID/Structures/Structures_test"
path_exp_file = "/home/malab/Desktop/BioEn-master/examples/scattering/files/experimental_data"
theta = "50" #theta
gl = "" #Optional- "" for loop

#Parameters for iBME
#will use exp_file
calc_rows_path = "{}/GP{}/calc_rows.txt" #path to calc_rows.txt (will be GP(n))
#will use theta
out_name = "{}/GP{}/"

##Create grid

#assign dro and r0 values and steps
#insert start value, end value + step, and step
dro_grid = np.arange(28, 30 + 2.0, 2.0) 
r0_grid =  np.arange(1.4, 1.6 + 0.1, 0.1)

grid_data = []
index = 0
for d in dro_grid:
    dro_val = round(d, 2)
    for r in r0_grid:
        r0_val = round(r, 2)
        grid_data.append([index, dro_val, r0_val])
        index += 1

grid_name = np.array(grid_data)
np.savetxt(
    "{}/grid_run.txt".format(sub_path), grid_name, fmt=["%d", "%.2f", "%.2f"],
    header="# d_rho r0", comments='')       

#Create dataframe of file names for later use
contents = pd.DataFrame(natsorted(os.listdir(path_structures)))

##Run do_gp

#move to new directory
env = os.environ.copy()
env["OUTPUT_DIR"] = sub_path

#structure path, experiment path, theta, gl (optional), grid document
run = subprocess.run(["../do_gp_v3.sh", path_structures, "{}/SASDLU4.dat".format(path_exp_file),
                      theta, gl, os.path.join(sub_path, "grid_run.txt"), sub_path], 
                      stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
print(f"STDOUT:\n {run.stdout}")
print(f"STDERR:\n {run.stderr}")

if run.returncode != 0:
    print(f"Script failed with return code {run.returncode}")
if run.returncode != 1:
    print(f"Script successful with return code {run.returncode}")

##File Modification

#Move GP files to correct subfolder location
#GP = range(0, 10000)
#folder = 0
#for i in GP:
#    if os.path.isdir(f"{main_path}/GP{i}"):
#        try:
#            shutil.move("{}/GP{}".format(main_path, i), sub_path)
#        except shutil.Error as e:
#            print(f"Error: {e}")
#    else:
#        print("All folders moved to sub path")
#        break
    
#Counter file loop and make columns
counter = range(0, 10000) #Arbitrary end- increase if more directories necessary
files = 0
for i in counter:
    if os.path.isdir("{}/GP{}".format(sub_path, i)):
        files += 1
        # Read calc_saxs.txt with frame index + intensities
        calc_rows = pd.read_csv("{}/GP{}/calc_saxs.txt".format(sub_path, i),
                             header=None, delim_whitespace=True)
     
        # Drop the first column (frame index)
        calc_rows = calc_rows.drop(columns=[0])

        # Write calc_rows.txt with header line # DATA=SAXS
        calc_path = '{}/GP{}/calc_rows.txt'.format(sub_path, i)
        with open(calc_path, 'w') as f:
            f.write("# DATA=SAXS\n")
        calc_rows.to_csv(calc_path, mode='a', header=False, index=False, sep=' ')
          
    else:
        print(f"Amount of GP directories:")
        break
print(files)

#Truncate File(s)
sim_length = len(calc_rows.columns)
print(sim_length)
exp_pd = pd.DataFrame(pd.read_csv("{}/SASDLU4.dat".format(path_exp_file),
                                  header=None, delim_whitespace=True))
exp_trun = exp_pd.iloc[:sim_length-1]
len_exp = len(exp_trun)
print(len_exp)
exp_trun.to_csv("{}/SASDLU4_trun.dat".format(path_exp_file),
                header=False, index=False, sep=' ')
with open(f"{path_exp_file}/SASDLU4_trun.dat", "r+") as f:
    content = f.read()
    f.seek(0, 0)
    f.write("# DATA=SAXS BOUNDS=UPPER\n" + content)

files = int(files)

##Grid to Dataframe for scan analysis
GRID_DF = pd.DataFrame(pd.read_csv("{}/grid_run.txt".format(sub_path), delim_whitespace=True))

#edit for solo grid run

##iBME 

for i in range(len(GRID_DF)):
    # Collect Parameters
    drho = GRID_DF.iloc[i]['d_rho']
    r0 = GRID_DF.iloc[i]['r0']

    #set parameters to nan initally
    chi2b = chi2a = phi = np.nan
    try:
        #Run iBME
        #truncated exp file, calc_rows, theta, output format
        run2 = iBME_script.iBMEf("{}/SASDLU4_trun.dat".format(path_exp_file),
                                 calc_rows_path.format(sub_path, i),
                                 theta, out_name.format(sub_path, i))
        

        
        logs = glob.glob(os.path.join(out_name.format(sub_path, i), "_ibme_*.log"))

        # Extract the number from the filename and sort numerically
        logs_sorted = sorted(logs, key=lambda x: int(re.search(r"_ibme_(\d+)\.log", x).group(1)))

        # Pick the last (largest iteration)
        log_file = logs_sorted[-1] if logs_sorted else None


        if log_file:
            with open(log_file) as lf:
                for L in lf:
                    if "CHI2 before optimization:" in L:   chi2b = float(L.split()[-1])
                    elif "CHI2 after optimization:" in L:  chi2a = float(L.split()[-1])
                    elif "Fraction of effective frames:" in L: phi = float(L.split()[-1])
        print(f"iBMEf run completed for GP{i}")
    except Exception as e:
        print(f"iBMEf failed for GP{i}: {e}")

    #paste collected parameters- nan will stay for failed iBME runs
    idx = (i)
    rows = []
    rows.append([idx, drho, r0, chi2b, chi2a, phi])

    #Create file with collected parameters
    grid = np.array(rows, dtype=float)
    np.savetxt("{}/GP{}/GRID_opt_{}".format(sub_path, i, i), grid,
       header="idx d_rho r0 CHI2_before CHI2_after PHI_eff", fmt="%.6g")
        


print("Success")

# Create summarized collected parameters file for plotting analysis

results = []
for i in range(len(GRID_DF)):
    folder_path = (f"{sub_path}/GP{i}/GRID_opt_{i}")
    if not os.path.exists(folder_path):
        print(f"Folder not found: {folder_path}")
        continue
    frames = pd.read_csv(folder_path)
    results.append(frames)
points = pd.concat(results, ignore_index=True)
print(points)
points.to_csv("{}/GRID_sum.txt".format(sub_path), index=False)
grid_file = "{}/GRID_sum.txt".format(sub_path)
  
##Plotting

# 1) load your grid file (change the filename if needed)
grid = np.loadtxt(grid_file)   

# 2) extract axes
dro = np.unique(grid[:,1])
r0  = np.unique(grid[:,2])

# 3) ensure rows are ordered (d_rho major, r0 minor) for reshape
order = np.lexsort((grid[:,1], grid[:,2]))  # sort by (dro, r0)
grid = grid[order]

# 4) ensure gamma exists (last column = ln(chi2_after / phi_eff))
chi2 = np.clip(grid[:,4], 1e-12, None)
phi  = np.clip(grid[:,5], 1e-12, None)
if grid.shape[1] < 7:
    gamma = np.log(chi2 / phi)
    grid = np.concatenate([grid, gamma[:,None]], axis=1)

# grid columns: [0]=idx [1]=d_rho [2]=r0 [3]=chi2_before [4]=chi2_after [5]=phi_eff [6]=gamma
chi2 = np.clip(grid[:,4], 1e-12, None)
phi  = np.clip(grid[:,5], 1e-12, None)
gam  = grid[:,6]  # already ln(chi2_after / phi_eff)

# reshape to (len(r0), len(dro)) for imshow with r0 on Y, d_rho on X
chi2_mat = np.log(chi2).reshape(len(r0), len(dro))
phi_mat  = (phi).reshape(len(r0), len(dro))
gam_mat  = gam.reshape(len(r0), len(dro))

# locate minimum gamma to highlight
min_y, min_x = np.unravel_index(np.nanargmin(gam_mat), gam_mat.shape)
best_dro = dro[min_x]
best_r0  = r0[min_y]

fig, axs = plt.subplots(1, 3, figsize=(18, 5), dpi=150)

im0 = axs[0].imshow(chi2_mat, origin='upper', aspect='auto')
axs[0].set_title(r'$\ln(\chi^2_{\mathrm{after}})$')
axs[0].scatter(min_x, min_y, s=60, marker='o', facecolors='none', edgecolors='k')
plt.colorbar(im0, ax=axs[0], fraction=0.046, pad=0.04,)

im1 = axs[1].imshow(phi_mat, origin='upper', aspect='auto')
axs[1].set_title(r'$\phi_{\mathrm{eff}}$')
axs[1].scatter(min_x, min_y, s=60, marker='o', facecolors='none', edgecolors='k')
plt.colorbar(im1, ax=axs[1], fraction=0.046, pad=0.04)

im2 = axs[2].imshow(gam_mat, origin='upper', aspect='auto')
axs[2].set_title(r'$\gamma=\ln(\chi^2_{\mathrm{after}}/\phi_{\mathrm{eff}})$')
axs[2].scatter(min_x, min_y, s=60, marker='o', facecolors='none', edgecolors='k')
plt.colorbar(im2, ax=axs[2], fraction=0.046, pad=0.04)

# tick labels: show every 2nd tick for readability
xticks = np.arange(0, len(dro), 2)
yticks = np.arange(0, len(r0), 2)
for ax in axs:
    ax.set_xticks(xticks); ax.set_xticklabels([f'{dro[i]:.2f}' for i in xticks], rotation=300)
    ax.set_yticks(yticks); ax.set_yticklabels([f'{r0[i]:.3f}' for i in yticks])
    ax.set_xlabel(r'$\delta\rho$  [$e/\mathrm{nm}^3$]')
axs[0].set_ylabel(r'$r_0/r_m$')

fig.suptitle(f'Best: δρ={best_dro:.2f}, r0={best_r0:.3f}', y=1.0)
plt.tight_layout()

## Save final summary of run

#Create Directory in /run_summary for each run
os.mkdir(f'{sub_path}/run_summary')
run_sum_path = f'{sub_path}/run_summary'
dir_num = 0
run_directory = "run_info_{}_{}"
info_path = os.path.join(run_sum_path, run_directory)
while os.path.isdir(info_path.format(today, dir_num)):
    dir_num += 1
info_path = info_path.format(today, dir_num)
os.mkdir(info_path)

#Amount of structures used
struc_counter = range(0, 10000) #Arbitrary end- increase if more directories necessary
struc = 0
for i in struc_counter:
    if os.path.isfile(f"{path_structures}/mm016_{i:03}.pdb"):
        struc += 1
    else:
        break

#Other parameters (can add more)
min_dro = min(grid[:,1])
max_dro = max(grid[:,1])
min_r0 = min(grid[:,2])
max_r0 = max(grid[:,2])

#Create file
summary = []
summary.append(f'Time of run: {cd}, Structure count: {struc}, Grid point count: {files}, Best dro: {best_dro}, Best r0: {best_r0}, Min dro: {min_dro}, Max dro: {max_dro}, Min r0: {min_r0}, Max r0: {max_r0}')
dfsummary = pd.DataFrame(summary)

#Modify name if required and save
run_num = 0
sumfile = os.path.join(info_path, "run_{}__{}.txt")
while os.path.isfile(sumfile.format(today, run_num)):
    run_num += 1
sumfile = sumfile.format(today, run_num)    
dfsummary.to_csv(sumfile, index=False)

heatmap = f'{info_path}/grid_heatmaps_{today}_{dir_num}.png'
fig.savefig(heatmap, dpi=300)
plt.show()

#Same for all grid points
pointsdf = pd.DataFrame(points)
gridrun_num = 0
gridsumfile = os.path.join(info_path, "gridrun_{}__{}.txt")
while os.path.isfile(gridsumfile.format(today, gridrun_num)):
    gridrun_num += 1 
gridsumfile = gridsumfile.format(today, gridrun_num)
pointsdf.to_csv(gridsumfile, index=False)

#Save pdb weights of most optimal grid point
weight_idx = GRID_DF.index[(GRID_DF['d_rho']== best_dro) & (GRID_DF['r0']== best_r0)].tolist()
weight_str = str(weight_idx[0])

copy_dir = f'{info_path}/GP{weight_str}'
shutil.copytree(f'{sub_path}/GP{weight_str}', copy_dir, dirs_exist_ok=True)

opt_weight = pd.DataFrame(pd.read_csv("{}/GP{}/_19.weights.dat".format(sub_path, weight_str), delim_whitespace=True, header=None))
opt_weight['2'] = opt_weight.iloc[:, 0].map(contents.iloc[:, 0])
opt_sorted = opt_weight.sort_values(by=1, ascending=False)
opt_sorted.to_csv(f'{info_path}/structure_weights_sorted_{today}_{dir_num}.txt', index=None)
