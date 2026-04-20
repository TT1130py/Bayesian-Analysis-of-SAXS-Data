import os
import glob
import subprocess
import pandas as pd
import numpy as np
import re
import shutil
import matplotlib.pyplot as plt
from datetime import date
from natsort import natsorted
from concurrent.futures import ProcessPoolExecutor, as_completed
import iBME_mod_SSH

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

#########--------- ARGUMENTS
#parser = argparse.ArgumentParser(description="Run iBME reweighting on fractions")
#parser.add_argument("save_path", type=str)
#parser.add_argument("dro_min", type=str)
#parser.add_argument("dro_max", type=str)
#parser.add_argument("r0_min", type=str)
#parser.add_argument("r0_max", type=str)

exp_path = "/home/tjaglal/Structures/Experiment/SASDLU4.dat"
struc_path = "/home/tjaglal/Structures/"
#args = parser.parse_args()
trun_path = "/home/tjaglal/Structures/Experiment/SASDLU4_trun.dat"
out_path = "/home/tjaglal/iBME/theta_scan"
today = date.today()
save_path = "/home/tjaglal/theta_scan_files"


######## PARAMS
dro_min = 38
dro_max = 38
r0_min = 1.34
r0_max = 1.36

########---------- Assign theta scan
theta_vals = np.array([10, 50, 100])
#########---------- FUNCTIONS

def ibme_worker(i, dro, r0, theta, calc_path, gp_out_dir, trun_path):
    os.makedirs(gp_out_dir, exist_ok=True)
    chi2b = chi2a = phi = np.nan

    try:
        # Run iBME
        iBME_mod_SSH.iBMEf(trun_path, calc_path, theta, f"{gp_out_dir}/")

        # Parse Logs
        logs = glob.glob(os.path.join(gp_out_dir, "_ibme_*.log"))
        logs_sorted = sorted(logs, key=lambda x: int(re.search(r"_ibme_(\d+)\.log", x).group(1)))
        log_file = logs_sorted[-1] if logs_sorted else None

        if log_file:
            with open(log_file) as lf:
                for L in lf:
                    if "CHI2 before optimization:" in L:
                        chi2b = float(L.split()[-1])
                    elif "CHI2 after optimization:" in L:
                        chi2a = float(L.split()[-1])
                    elif "Fraction of effective frames:" in L:
                        phi = float(L.split()[-1])
        print(f"GP{i} optimized.")
    except Exception as e:
        print(f"iBME failed for GP{i}: {e}")

    rows = [[i, dro, r0, chi2b, chi2a, phi]]
    grid = np.array(rows, dtype=float)
    np.savetxt(os.path.join(gp_out_dir, f"GRID_opt_{i}"), grid, header="idx d_rho r0 CHI2_before CHI2_after PHI_eff",
               fmt="%.6g")

    return {"idx": i, 'd_rho': dro, "r0": r0, "CHI2_before": chi2b, "CHI2_after": chi2a, "PHI": phi}

#######------- MAIN

#Concat fractions in main structure folder
if __name__ == "__main__":
    grid_file_path = os.path.join(save_path, "grid_full.txt")
    GRID_DF = pd.read_csv(grid_file_path, sep='\s+', header=None, names=['index', 'dro', 'r0'])
    
    sample_saxs = os.path.join(save_path, "GP0_all_saxs.txt")
    with open(sample_saxs, 'r') as f:
        sim_length = len(f.readline().split())
    
    exp_pd = pd.read_csv(exp_path, header=None, sep='\s+')
    exp_trun = exp_pd.iloc[:sim_length-1]
    exp_trun.to_csv(trun_path, header=False, index=False, sep=' ')
    
    with open(trun_path, "r+") as f:
        content = f.read()
        f.seek(0, 0)
        f.write("# DATA=SAXS BOUNDS=UPPER\n" + content)
    
    all_chi2 = []
    all_skl = []
    print(f"Starting iBME optimizations for {len(GRID_DF)} grid points...")
    
    for theta in theta_vals:
        run_fol = f"iBME_dro_{dro_min}_to_{dro_max}_r0_{r0_min}_to_{r0_max}_theta_{theta}"
        ibme_out_dir = os.path.join(out_path, run_fol)
        os.makedirs(ibme_out_dir, exist_ok=True)
    
        results = []
    
        with ProcessPoolExecutor() as executor:
            futures = []
            for i in range(len(GRID_DF)):
                dro = GRID_DF.iloc[i]['dro']
                r0 = GRID_DF.iloc[i]['r0']
                calc_path = os.path.join(save_path, f"GP{i}_all_saxs.txt")
                gp_out_dir = os.path.join(ibme_out_dir, f"GP{i}")
    
                #Paths here
                futures.append(executor.submit(ibme_worker, i, dro, r0, theta, calc_path, gp_out_dir, trun_path))
    
            for future in as_completed(futures):
                res = future.result()
                if 'error' not in res:
                    results.append(res)
                else:
                    print(f"iBME failed for GP{res['idx']}: {res['error']}")
    
        ######------- Scan theta folders to analyze CHI2 and phi/Skl for plotting
    
        points = pd.DataFrame(results)
        grid_sum_path = os.path.join(ibme_out_dir, "GRID_sum.txt")
        points.to_csv(grid_sum_path, index=False)
    
        grid = np.loadtxt(grid_sum_path, skiprows=1, delimiter=',')  # Skip pandas header
    
        # Extract dRho and r0 values from the grid
        dro_vals = np.unique(grid[:, 1])
        r0_vals = np.unique(grid[:, 2])
        order = np.lexsort((grid[:, 1], grid[:, 2]))
        grid = grid[order]
    
        # Extract CHI2 and phi values
        chi2 = np.clip(grid[:, 4], 1e-12, None)
        phi = np.clip(grid[:, 5], 1e-12, None)
    
        # Convert all phi values to Skl
        skl = -np.log(phi)
    
        # Find gamma using formula transformation
        gamma = np.log(chi2) + skl
    
        # reshape for minimum values
        chi2_mat = np.log(chi2).reshape(len(r0_vals), len(dro_vals))
        skl_mat = skl.reshape(len(r0_vals), len(dro_vals))
        gam_mat = gamma.reshape(len(r0_vals), len(dro_vals))
    
        # find the best SAXS parameters
        min_y, min_x = np.unravel_index(np.nanargmin(gam_mat), gam_mat.shape)
        best_dro = dro_vals[min_x]
        best_r0 = r0_vals[min_y]
    
    
        all_chi2.append(chi2_mat[min_y, min_x])
        all_skl.append(skl_mat[min_y, min_x])
    
    #####-------- Plotting
    chi_np = np.array(all_chi2)
    skl_np = np.array(all_skl)
    
    #DataFrame for organized plotting
    all_data = pd.DataFrame({"Theta": theta_vals, "chi": chi_np, "skl": skl_np})
    
    #Plot
    fig = plt.figure(figsize=[12, 8])
    ax = fig.add_subplot(111)
    
    for theta in range(len(theta_vals)):
        chi_plot = chi_np[theta]
        skl_plot = skl_np[theta]
        ax.scatter(chi_plot, skl_plot, marker='o', label=theta_vals[theta])
    
    plt.grid()
    ax.set_xlabel("chi2")
    ax.set_ylabel("skl")
    ax.legend(ncol=2)
    
    fig_path = os.path.join(out_path, "L_curve.png")
    fig.savefig(fig_path, dpi=300)


