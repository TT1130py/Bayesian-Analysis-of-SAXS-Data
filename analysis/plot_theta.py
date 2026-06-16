import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Match the paths and params from your original run
out_path = "/users/t/j/tjaglal/Projects/iBME/theta_scan/30_50_1.44_1.62_1000_10000"
dro_min = 30
dro_max = 50
r0_min = 1.44
r0_max = 1.64
theta_vals = np.array([1000, 2000, 4000, 6000, 8000, 10000])

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
