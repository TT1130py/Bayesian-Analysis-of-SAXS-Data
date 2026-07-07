#####----- Examine fitting based on PEPSI simulations
import numpy as np
import argparse
import os
import glob
from natsort import natsorted
import re
import yaml
from pathlib import Path

#####----- FUNCTIONS

def extract_chi2(file_path):
    chi2_values = []
    chi_pattern = re.compile(r"Chi\^2\.+\s*:\s*([\d\.e\+\-]+)")

    with open(file_path, 'r') as file:
        for line in file:
            match = chi_pattern.search(line)
            if match:
                chi2_value = float(match.group(1))
                chi2_values.append(chi2_value)

    return chi2_values

def iterate_dir(type, sim_path):
    chi_avg = []
    for gp in gps:
        if type == "c_term":
            sim_dirs = glob.glob(os.path.join(sim_path, 'MD*_*'))
        else:
            sim_dirs = glob.glob(os.path.join(sim_path, 'mm*'))
        sorted_sim_dirs = natsorted(sim_dirs)

        gp_chi = []
        for data_file in sorted_sim_dirs:
            gp_dir = os.path.join(data_file, "GP{}".format(gp))
            log_file = os.path.join(gp_dir, "logPEPSI")
            chi2_ind = extract_chi2(log_file)
            gp_chi.append(chi2_ind)
        mean_chi2 = np.mean(gp_chi)
        chi_avg.append(mean_chi2)
    all_chi = np.array(chi_avg)

    return all_chi

def best_fit(all_chi):
    min_index = np.argmin(all_chi)
    min_chi = all_chi[min_index]

    print(f"The best fit grid point is GP = {min_index} with an average $\chi^2$ = {min_chi}")

    return min_index

def structure_fit(min_gp, type):
    results = []
    current_filename = None

    pdb_pattern = re.compile(r"PDB filename\.+\s*:\s*(.+)")
    chi2_pattern = re.compile(r"Chi\^2\.+\s*:\s*([\d\.e\+\-]+)")

    if type == "c_term":
        sim_dirs = glob.glob(os.path.join(sim_path, 'MD*_*'))
    else:
        sim_dirs = glob.glob(os.path.join(sim_path, 'mm*'))
    sorted_sim_dirs = natsorted(sim_dirs)

    for data_file in sorted_sim_dirs:
        gp_dir = os.path.join(data_file, "GP{}".format(min_gp))
        log_file = os.path.join(gp_dir, "logPEPSI")
        with open(log_file, 'r') as file:
            for line in file:
                pdb_match = pdb_pattern.search(line)
                if pdb_match:
                    full_path = pdb_match.group(1).strip()
                    current_filename = os.path.basename(full_path)

                chi2_match = chi2_pattern.search(line)
                if chi2_match:
                    chi2_value = float(chi2_match.group(1))

                    if current_filename:
                        results.append({
                            "filename": current_filename,
                            "chi2": chi2_value
                        })
                        current_filename = None

    return results

def eval_structures(best_structures, cutoff, min_gp):
    top = best_structures[:cutoff]
    with open(summary_file, "w") as out_file:
        # Write a clean, human-readable table header
        out_file.write(f"--- Top {cutoff} Best Fitting Pepsi-SAXS Structures at grid point {min_gp} ---\n")
        out_file.write(
            f"{'Rank':<6}{'File Name':<25}{'Chi2':<10}\n"
        )
        out_file.write("-" * 65 + "\n")

        for rank, entry in enumerate(top, 1):
            out_file.write(
                f"{rank:<6}"
                f"{entry['filename']:<25}"
                f"{entry['chi2']:<10.2f}"
            )

def main():
    #####----- Initialize CLI arguments and configuration
    parser = argparse.ArgumentParser()

    parser.add_argument("--config", type=str, default="config.yaml", help="Path to main YAML file")
    args = parser.parse_parser_args() if hasattr(parser, 'parse_parser_args') else parser.parse_args()

    with open(args.config, "r") as f:
        master_config = yaml.safe_load(f)

    pi_config = master_config.get("pepsi_iteration", {})
    if not pi_config:
        print("pepsi_fit config not found")
        return

    sim_path = pi_config.get("sim_path", "")
    summary_file = pi_config.get("summary_file", "")
    gps = int(pi_config.get("gps", ""))
    cutoff = int(pi_config.get("cutoff", ""))
    type = pi_config.get("type", "")

    all_chi = iterate_dir(type, sim_path)
    min_gp = best_fit(all_chi)
    result = structure_fit(min_gp, type)
    eval_structures(result, cutoff, gps)