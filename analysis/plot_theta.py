import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import yaml

#####----- CLI Config and Initialization
parser= argparse.ArgumentParser()
parser.add_argument("--config", type=str, default="config.yaml", help="Path to main yaml file")
args = parser.parse_parser_args() if hasattr(parser, 'parse_parser_args') else parser.parse_args()

with open(args.config, "r") as f:
    master_config = yaml.safe_load(f)

pt_config = master_config.get("plot_theta", {})
if not pt_config:
    print("plot theta config not found")
    return

out_path = pt_config.get["out_path", ""]
dro_min = pt_config.get["dro_min", ""]
dro_max = pt_config.get["dro_max", ""]
r0_min = pt_config.get["r0_min", ""]
r0_max = pt_config.get["r0_max", ""]
theta_val_1 = pt_config.get["theta_val_1", ""]
theta_val_2 = pt_config.get["theta_val_2", ""]
theta_val_3 = pt_config.get["theta_val_3", ""]
theta_val_4 = pt_config.get["theta_val_4", ""]
theta_val_5 = pt_config.get["theta_val_5", ""]
theta_val_6 = pt_config.get["theta_val_6", ""]

#####----- MAIN

all_chi2 = []
all_skl = []
valid_thetas = [] # Keeps track of which thetas actually have a saved file

for theta in theta_vals:
    # Reconstruct the folder path where the data is saved
    run_fol = f"iBME_dro_{dro_min}_to_{dro_max}_r0_{r0_min}_to_{r0_max}_theta_{theta}"
    grid_sum_path = os.path.join(out_path, run_fol, "GRID_sum.txt")
    
    if os.path.exists(grid_sum_path):
        # 1. Read the saved dataframe
        df = pd.read_csv(grid_sum_path)
        
        # 2. Extract arrays and clip them to avoid log(0)
        chi2 = np.clip(df['CHI2_after'].values, 1e-12, None)
        phi = np.clip(df['PHI'].values, 1e-12, None)
        
        # 3. Calculate skl and gamma
        skl = -np.log(phi)
        gamma = np.log(chi2) + skl
        
        # 4. Find the row with the minimum gamma
        best_idx = np.nanargmin(gamma)
        
        # 5. Append the best chi2 and skl to our lists
        all_chi2.append(chi2[best_idx])
        all_skl.append(skl[best_idx])
        valid_thetas.append(theta)
        
        # Optional: Print to terminal so you can see the exact values
        print(f"Theta: {theta:3} | Best Chi2: {chi2[best_idx]:.4f} | Best Skl: {skl[best_idx]:.4f}")
    else:
        print(f"File not found, skipping: {grid_sum_path}")

# --- Plotting ---
if valid_thetas:
    fig = plt.figure(figsize=[12, 8])
    ax = fig.add_subplot(111)
    
    for i in range(len(valid_thetas)):
        ax.scatter(all_skl[i], all_chi2[i], marker='o', s=300, label=f"Theta: {valid_thetas[i]}")
    
    # Optional: Draw a line connecting the points to make the L-curve easier to see
    ax.plot(all_skl, all_chi2, linestyle='-', color='blue', alpha=0.5)

    plt.grid()
    ax.set_xlabel("skl")
    ax.set_ylabel("chi2")
    ax.set_title("Skl vs Chi^2 for Optimzed Thetas")
    ax.legend(ncol=2)
    
    fig_path = os.path.join(out_path, "L_curve_replot.png")
    fig.savefig(fig_path, dpi=300)
    print(f"\nSuccess! Plot saved to {fig_path}")
else:
    print("No data found to plot. Check your directory paths.")
