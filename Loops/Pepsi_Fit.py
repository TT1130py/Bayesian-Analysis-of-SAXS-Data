#####----- PEPSI_FIT
import os
import glob
import yaml
import argparse
import subprocess
from pathlib import Path
from natsort import natsorted
import re

#####----- REPLACE WITH ARGPARSE HERE ONCE FINALIZED
structure_path = ""
output_fit = ""
path_exp_file = ""
path_to_pepsi = ""
summary_file = ""
cutoff = ""

#####----- FUNCTIONS
def pepsi_loop(structure_path, path_to_pepsi, path_exp_file):
    counter = 0
    unordered_files = [f for f in structure_path.iterdir() if f.is_file()]
    file_list = natsorted(unordered_files, key=lambda p: p.name)
    for i in file_list:
        result = subprocess.run(["{}/Pepsi-SAXS".format(path_to_pepsi), os.path.join(structure_path, i.name), path_exp_file],
                                 capture_output=True, text=True)

        #Counter
        counter += 1

        # Error
        if result.returncode != 0:  # error
            print(f" pepsi failed for {result}")
            print(result.stderr)
            continue

    print(f"{counter} pepsi loops completed")

def find_fits():
    fit_extension = ("*.fit")

    metrics_pattern = re.compile(
        r"r0:\s+([0-9.-]+)\s+d_rho:\s+([0-9.-]+)\s+Chi2:\s+([0-9.-]+)")

    best_structures = []

    unordered_files = [f for f in structure_path.glob(fit_extension) if f.is_file()]
    file_list = natsorted(unordered_files, key=lambda p: p.name)
    for i in file_list:
        try:
            with open(i, "r") as f:
                for _ in range(10):
                    line = f.readline()
                    if not line:
                        break

                    match = metrics_pattern.search(line)
                    if match:
                        r0_val = float(match.group(1))
                        d_rho_val = float(match.group(2))
                        chi2_val = float(match.group(3))

                        # Store the filename and its corresponding chi2 score
                        best_structures.append(
                            {"file": i.name, "chi2": chi2_val, "d_rho": d_rho_val, "r0": r0_val}
                        )
                        break
        except Exception as e:
            print(f"Error reading {i.name}: {e}")

    best_structures.sort(key=lambda x: x["chi2"])
    return best_structures

def rank_fits(best_structures, cutoff):
    top = best_structures[:cutoff]
    with open(summary_file, "w") as out_file:
        # Write a clean, human-readable table header
        out_file.write(f"--- Top {cutoff} Best Fitting Pepsi-SAXS Structures ---\n")
        out_file.write(
            f"{'Rank':<6}{'File Name':<25}{'Chi2':<10}{'r0':<10}{'d_rho':<10}\n"
        )
        out_file.write("-" * 65 + "\n")

        for rank, entry in enumerate(top, 1):
            out_file.write(
                f"{rank:<6}"
                f"{entry['file']:<25}"
                f"{entry['chi2']:<10.2f}"
                f"{entry['r0']:<10.2f}"
                f"{entry['d_rho']:<10.2f}\n"
            )

def main():
    pepsi_loop(structure_path, path_to_pepsi, path_exp_file)
    sorted_structures = find_fits(structure_path)
    rank_fits(sorted_structures, cutoff)

if __name__ == "__main__":
    main()