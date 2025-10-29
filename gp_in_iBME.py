#gp in iBME
#running do_gp.sh from inside a python script rather than iBME in shell script
#RUN FROM VENV

import subprocess
import os
import iBME_script
import pandas as pd
import glob 
import numpy as np
import re
import matplotlib.pyplot as plt


##Set Parameters

#Working directory
print("Current working directory: {0}".format(os.getcwd()))
os.chdir('/home/malab/iBME')
print("Current working directory: {0}".format(os.getcwd()))

#Parameters for do_gp
path_structures = "/home/malab/Desktop/PDB_ID/Structures/Structures_test"
path_exp_file = "/home/malab/Desktop/BioEn-master/examples/scattering/files/experimental_data"
theta = "50" #theta
gl = "" #Optional- "" for loop
path_gpdoc = "/home/malab/iBME" #path to grid 

#Parameters for iBME
#will use exp_file
calc_rows_path = "/home/malab/iBME/GP{}/calc_rows.txt" #path to calc_rows.txt (will be GP(n))
#will use theta
out_name = "/home/malab/iBME/GP{}/"

##Run do_gp

#structure path, experiment path, theta, gl (optional), grid document
run = subprocess.run(["./do_gp_v3.sh", path_structures, "{}/SASDLU4.dat".format(path_exp_file),
                      theta, gl, "{}/GRID_tau".format(path_gpdoc)], 
                      stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
print(f"STDOUT:\n {run.stdout}")
print(f"STDERR:\n {run.stderr}")

if run.returncode != 0:
    print(f"Script failed with return code {run.returncode}")
if run.returncode != 1:
    print(f"Script successful with return code {run.returncode}")

##File Modification

#Counter file loop and make columns
counter = range(0, 10000) #Arbitrary end- increase if more directories necessary
files = 0
for i in counter:
    if os.path.isdir("/home/malab/iBME/GP{}".format(i)):
        files += 1
        # Read calc_saxs.txt with frame index + intensities
        calc_rows = pd.read_csv("/home/malab/iBME/GP{}/calc_saxs.txt".format(i),
                             header=None, delim_whitespace=True)
     
        # Drop the first column (frame index)
        calc_rows = calc_rows.drop(columns=[0])

        # Write calc_rows.txt with header line # DATA=SAXS
        calc_path = '/home/malab/iBME/GP{}/calc_rows.txt'.format(i)
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
GRID_DF = pd.DataFrame(pd.read_csv("{}/GRID_tau".format(path_gpdoc), delim_whitespace=True))

##iBME 


for i in range(files):
    try:
        #Run iBME
        #truncated exp file, calc_rows, theta, output format
        out_dir = out_name.format(i)
        run2 = iBME_script.iBMEf("{}/SASDLU4_trun.dat".format(path_exp_file), 
                                 calc_rows_path.format(i),
                                 theta, out_name.format(i))
        
        #Collect Parameters 
        drho = GRID_DF.iloc[i]['d_rho']
        r0 = GRID_DF.iloc[i]['r0']
        
        logs = glob.glob(os.path.join(out_name.format(i), "_ibme_*.log"))

        # Extract the number from the filename and sort numerically
        logs_sorted = sorted(logs, key=lambda x: int(re.search(r"_ibme_(\d+)\.log", x).group(1)))

        # Pick the last (largest iteration)
        log_file = logs_sorted[-1] if logs_sorted else None

        chi2b = chi2a = phi = np.inf
        if log_file:
            with open(log_file) as lf:
                for L in lf:
                    if "CHI2 before optimization:" in L:   chi2b = float(L.split()[-1])
                    elif "CHI2 after optimization:" in L:  chi2a = float(L.split()[-1])
                    elif "Fraction of effective frames:" in L: phi = float(L.split()[-1])
        idx = (i)
        rows = []
        rows.append([idx, drho, r0, chi2b, chi2a, phi])

        #Create file with collected parameters
        grid = np.array(rows, dtype=float)
        np.savetxt("/home/malab/iBME/GP{}/GRID_opt_{}".format(i, i), grid,
           header="idx d_rho r0 CHI2_before CHI2_after PHI_eff", fmt="%.6g")
        
        print(f"iBMEf run completed for GP{i}")
    except Exception as e:
       print(f"iBMEf failed for GP{i}: {e}")

print("Success")

# Create summarized collected parameters file for plotting analysis

results = []
for i in range(files):
    frames = pd.DataFrame(pd.read_csv(f"/home/malab/iBME/GP{i}/GRID_opt_{i}"))
    results.append(frames)
points = pd.concat(results)
print(points)
points.to_csv("/home/malab/iBME/GRID_sum.txt", index=False)
grid_file = "/home/malab/iBME/GRID_sum.txt"
  
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
plt.colorbar(im0, ax=axs[0], fraction=0.046, pad=0.04)

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

fig.suptitle(f'Best: δρ={best_dro:.2f}, r0={best_r0:.3f}', y=1.02)
plt.tight_layout()
# plt.savefig('grid_heatmaps.png', dpi=300)
plt.show()
