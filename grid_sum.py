import os
import glob
import argparse
import pandas as pd
import numpy as np
import re


def recover_grid_sum():
    parser = argparse.ArgumentParser(description="Recover GRID_sum.txt from partial iBME runs")
    parser.add_argument("save_path", type=str, help="Path containing grid_full.txt and iBME_results folder")
    args = parser.parse_args()

    grid_file_path = os.path.join(args.save_path, "grid_full.txt")
    ibme_out_dir = os.path.join(args.save_path, "iBME_results")

    try:
        GRID_DF = pd.read_csv(grid_file_path, sep='\s+', header=None, names=['index', 'dro', 'r0'])
    except FileNotFoundError:
        print(f"Error: Could not find master grid at {grid_file_path}")
        return

    results = []
    missing_count = 0

    print(f"Scanning {len(GRID_DF)} grid points...")

    for _, row in GRID_DF.iterrows():
        idx = int(row['index'])
        dro = row['dro']
        r0 = row['r0']

        gp_out_dir = os.path.join(ibme_out_dir, f"GP{idx}")

        #Initialize metrics as NaN
        chi2b = np.nan
        chi2a = np.nan
        phi = np.nan

        #Look for logs if the GP folder exists
        if os.path.isdir(gp_out_dir):
            logs = glob.glob(os.path.join(gp_out_dir, "_ibme_*.log"))

            if logs:
                #Sort logs to pull data from the most recent run
                logs_sorted = sorted(logs, key=lambda x: int(re.search(r"_ibme_(\d+)\.log", x).group(1)) if re.search(
                    r"_ibme_(\d+)\.log", x) else 0)
                log_file = logs_sorted[-1]

                try:
                    with open(log_file, 'r') as lf:
                        for L in lf:
                            if "CHI2 before optimization:" in L:
                                chi2b = float(L.split()[-1])
                            elif "CHI2 after optimization:" in L:
                                chi2a = float(L.split()[-1])
                            elif "Fraction of effective frames:" in L:
                                phi = float(L.split()[-1])
                except Exception as e:
                    print(f"Warning: Could not read log for GP{idx}: {e}")
                    missing_count += 1
            else:
                missing_count += 1
        else:
            missing_count += 1

        #Append the row, preserving NaN where logs were missing
        results.append([idx, dro, r0, chi2b, chi2a, phi])

    #Reconstruct the final dataframe
    recovered_df = pd.DataFrame(results, columns=["idx", "d_rho", "r0", "CHI2_before", "CHI2_after", "PHI_eff"])

    #Save to CSV in the exact format of the original script
    grid_sum_path = os.path.join(ibme_out_dir, "GRID_sum.txt")
    recovered_df.to_csv(grid_sum_path, index=False)

    print(f"\nRecovery Complete!")
    print(f"Saved to: {grid_sum_path}")
    print(f"Total points: {len(recovered_df)}")
    print(f"Points missing iBME output (NaN applied): {missing_count}")


if __name__ == "__main__":
    recover_grid_sum()