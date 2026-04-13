from os import WCONTINUED

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from natsort import natsort_keygen
import argparse
import os

system = "VACC" #Either VACC or Local
##########------ PLOT SAXS CURVE COMPARISON OF EXPERIMENT VS WEIGHTED AVERAGE SIMULATION

##########------ ARGUMENTS AND ABSOLUTE PATHS
if system == "VACC":
    print("VACC")
    #parser = argparse.ArgumentParser(description='Fit Calculated SAXS Ensemble to experiment')
    #parser.add_argument("dro", type=float)
    #parser.add_argument("r0", type=float)
    #arser.add_argument("save_path", type=str)
    #args = parser.parse_args()
elif system == "Local":
    print("LOCAL")
    #sim_path = "{}/iBME_results/structure_weights_sorted_*.txt".format(args.save_path)
else:
    print("System not supported")

#path to the experimental SAXS data
path_exp_file = "/home/malab/Desktop/BioEn-master/examples/scattering/files/experimental_data"
#path to the "structure weights sorted" file that contains structure file name and its weight (output from iBME)

#path to the simulated SAXS curves of above files
structure_path = "/home/malab/iBME/pepsi_output"

##########------ FUNCTIONS
def match_files(sim_file):
    sim_pd = pd.read_csv(sim_file, sep='\\t', header=0)
    sim_pd["PDB_Name"] = sim_pd["PDB_Name"].str.replace('.pdb','',regex=False)

    grid_df = pd.read_csv("/home/malab/Downloads/curve_test/grid_full.txt", sep='\s+', header=None, names=['index', 'dro', 'r0'])
    gp = grid_df.loc[(grid_df['dro'] == 30) & (grid_df['r0'] == 1.22), 'index'].iloc[0]
    #grid_df = pd.read_csv(os.path.join(args.save_path, "grid_full.txt"), sep='\s+', header=None, names=['index', 'dro', 'r0'])
    #gp = grid_df.loc[(grid_df['dro'] == args.dro) & (grid_df['r0'] == args.r0), 'index'].iloc[0]

    #compiled_saxs = np.genfromtxt("{}/compiled_GPs/GP{}_all_saxs.txt".format(args.save_path, str(gp)))
    compiled_saxs = np.genfromtxt("/home/malab/Downloads/curve_test/GP28_all_saxs.txt")
    compiled_df = pd.DataFrame(compiled_saxs).reset_index(drop=True)
    iq_sim_matrix = pd.DataFrame(compiled_df.drop(columns=[0]).values)
    pdb_sorted = sim_pd.sort_values(by="PDB_Name", key=natsort_keygen()).reset_index(drop=True)

    pd.set_option('display.float_format', '{:.4e}'.format)

    sorted_weights = pdb_sorted.iloc[:,1]
    s_full = np.genfromtxt("/home/malab/Downloads/curve_test/qvals.txt", skip_header=0)
    s_values = s_full[:iq_sim_matrix.shape[1]]
    print("Breakpt")
    return s_values, iq_sim_matrix, sorted_weights, pdb_sorted["PDB_Name"]

def VACC_average_curve(sim_file, pdb_names, s_val, iq_val):
    sim_pd = pd.read_csv(sim_file, sep='\\t', header=0)
    sim_pd["PDB_Name"] = sim_pd["PDB_Name"].str.replace('.pdb','',regex=False)

    weight_map = dict(zip(sim_pd['PDB_Name'], sim_pd['1']))
    ordered_weights = np.array([weight_map.get(name, 0.0) for name in pdb_names])
    iq_array = np.array(iq_val)

    weighted_matrix = iq_array * ordered_weights[:, np.newaxis]
    avg_iq = np.sum(weighted_matrix, axis=0)

    sim_merge = pd.concat([pd.DataFrame(s_val), pd.DataFrame(avg_iq)], axis=1)
    print("Breakpt")

    return len(sim_merge), sim_merge.iloc[:,0], sim_merge.iloc[:,1]
def simulated_curves(sim_file):
    sim_pd = pd.read_csv(sim_file, header=0)
    sim_pd["PDB_Name"] = sim_pd["PDB_Name"].str.replace('.pdb','',regex=False)

    all_weights = sim_pd[1]
    curves = len(sim_pd)
    all_sim_s = []
    all_sim_iq = []
    filenames = []

    for i in range(curves):
        if "mm" in sim_pd.iloc[i,2]:
            file = sim_pd.iloc[i,2]
        else:
            continue
        read_file = pd.read_csv("{}/{}.out".format(structure_path, file), sep=r"\s+", skiprows=5)
        s_sim = read_file.iloc[:,0] * 10
        iq_sim = read_file.iloc[:,1]

        all_sim_s.append(s_sim)
        all_sim_iq.append(iq_sim)
        filenames.append(file)

    df_s = pd.DataFrame(all_sim_s, index=filenames)
    df_iq = pd.DataFrame(all_sim_iq, index=filenames)

    fixed_weights = all_weights.drop(all_weights.index[0]).reset_index(drop=True)
    df_weights = pd.DataFrame({"weights": fixed_weights.values, "structure file": filenames})

    length = s_sim.iloc[-1]

    return length, df_s, df_iq, file, df_weights, all_weights

def average_curve(s_array, iq_array, weights):
    s_values = s_array.iloc[0, :]
    iq_tran = iq_array.transpose()
    merge = pd.concat([s_values, iq_tran], axis=1)

    merge.columns.values[0] = "s"
    weighted_merge = merge.copy()

    for i in range(len(merge.columns)):
        if "mm" in merge.columns[i]:
            structure = merge.columns[i]
        else:
            continue

        match = weights.loc[weights["structure file"] == structure, "weights"]
        if match.empty:
            w = 0
        else:
            w = match.iloc[0]
            if pd.isna(w):
                w = 0

        for x in range(len(merge)):
            iq_i = merge.loc[x, structure]
            weighted_merge.loc[x, structure] = iq_i * w

    averaged_merge = pd.DataFrame({"s": s_values, "iq": np.nan})
    for i in range(len(merge)):
        series = weighted_merge.iloc[i, 1:]
        series_np = np.array(series)
        average = np.sum(series_np)
        averaged_merge.loc[i, "iq"] = average

    return averaged_merge.loc[:,"s"], averaged_merge.loc[:,"iq"]

def experimental_curve(path_exp_file, sim_length):
    exp_pd = pd.DataFrame(pd.read_csv("{}/SASDLU4.dat".format(path_exp_file),
                                      header=None, sep=r"\s+"))

    s = exp_pd.loc[:,0]
    iq = exp_pd.loc[:,1]
    err = exp_pd.loc[:,2]

    s_trun = s.iloc[:sim_length].values
    iq_trun = iq.iloc[:sim_length].values
    err_trun = err.iloc[:sim_length].values

    plt.errorbar(s, iq, yerr=err, fmt='o', ecolor="lightgray", markersize=3)
    plt.yscale("log")
    plt.ylabel("i(q)")
    plt.xlabel("s")
    plt.title("Experimental SAXS curve with x log")
    plt.show()

    return s_trun, iq_trun, err_trun, s, iq, err
    print("Breakpt")

def plot_compare(s_sim, iq_sim, s, iq, err, s_full, iq_full, err_full, f_name):
    scale = np.sum(iq * iq_sim) / np.sum(iq_sim**2)

    scaled_iq_sim = iq_sim * scale

    fig, ax = plt.subplots(figsize = (10,10))
    ax.errorbar(s, iq, yerr = err, fmt= 'o', markersize=3, ecolor="lightgray", label="Experiment")
    ax.set_yscale("log")

    ax.plot(s_sim, scaled_iq_sim, zorder=3, lw=3, label=f_name)
    ax.set_ylabel("i(q)")
    ax.set_xlabel("s")
    ax.set_title("Simulated SAXS fit with Experiment - truncated")
    ax.legend()
    plt.show()

    fig_2, ax_2 = plt.subplots(figsize = (10,10))
    ax_2.errorbar(s_full, iq_full, yerr = err_full, fmt= 'o', markersize=3, ecolor="lightgray", label="Experiment")
    ax_2.set_yscale("log")

    ax_2.plot(s_sim, scaled_iq_sim, zorder=3, lw=3, label=f_name)
    ax_2.set_ylabel("i(q)")
    ax_2.set_xlabel("s")
    ax_2.set_title("Simulated SAXS fit with Experiment - full")
    ax.legend()
    plt.show()

##########------ MAIN

def main(run):
    if run == "VACC":
        #Match Files
        s, concat_merge, weights, f_name = match_files("/home/malab/Downloads/curve_test/structure_weights_sorted_38.0_1.34_50.0.csv")

        lent, s_weighted, iq_weighted = VACC_average_curve("/home/malab/Downloads/curve_test/structure_weights_sorted_38.0_1.34_50.0.csv", f_name, s, concat_merge)

        # Create the SAXS curves from experimental data for plotting
        angletrun, intensetrun, errtrun, anglefull, intensefull, errfull = experimental_curve(path_exp_file, lent)

        # Plot experiment vs weighted average simulation
        plot_compare(s_weighted, iq_weighted, angletrun, intensetrun, errtrun, anglefull, intensefull, errfull, "GP0_all_saxs")

    elif run == "Local":
        #Find all simulated SAXS curves and place in dataframe
        lent, sangle, intense, file_name, weights_df, wlist = simulated_curves(sim_path)

        # Perform the weighted average SAXS curve
        s_weighted, iq_weighted = average_curve(sangle, intense, weights_df)

        # Create the SAXS curves from experimental data for plotting
        angletrun, intensetrun, errtrun, anglefull, intensefull, errfull = experimental_curve(path_exp_file, lent)

        # Plot experiment vs weighted average simulation
        plot_compare(s_weighted, iq_weighted, angletrun, intensetrun, errtrun, anglefull, intensefull, errfull,
                     file_name)
    else:
        print("System name error")


main(system)
print("Breakpt")