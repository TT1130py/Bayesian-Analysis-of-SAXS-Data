import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

##########------ PLOT SAXS CURVE COMPARISON OF EXPERIMENT VS WEIGHTED AVERAGE SIMULATION

#path to the experimental SAXS data
path_exp_file = "/home/malab/Desktop/BioEn-master/examples/scattering/files/experimental_data"
#path to the "structure weights sorted" file that contains structure file name and its weight (output from iBME)
sim_path = "/home/malab/iBME/run_2025-11-20_0/run_summary/run_info_2025-11-20_0/structure_weights_sorted_2025-11-20_0.txt"
#path to the simulated SAXS curves of above files
structure_path = "/home/malab/iBME/pepsi_output"

##########------ FUNCTIONS

def simulated_curves(sim_file):
    sim_pd = pd.read_csv(sim_file, header=None)
    sim_pd[2] = sim_pd[2].str.replace('.pdb','',regex=False)

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

    mask = s < sim_length
    s_trun = s[mask]
    iq_trun = iq[mask]
    err_trun = err[mask]

    plt.errorbar(s, iq, yerr=err, fmt='o', ecolor="lightgray", markersize=3)
    plt.yscale("log")
    plt.ylabel("i(q)")
    plt.xlabel("s")
    plt.title("Experimental SAXS curve with x log")
    plt.show()

    return s_trun, iq_trun, err_trun, s, iq, err
    print("Breakpt")

def plot_compare(s_sim, iq_sim, s, iq, err, s_full, iq_full, err_full, f_name):
    fig, ax = plt.subplots(figsize = (10,10))
    ax.errorbar(s, iq, yerr = err, fmt= 'o', markersize=3, ecolor="lightgray", label="Experiment")
    ax.set_yscale("log")

    ax.plot(s_sim, iq_sim, zorder=3, lw=3, label=f_name)
    ax.set_ylabel("i(q)")
    ax.set_xlabel("s")
    ax.set_title("Simulated SAXS fit with Experiment - truncated")
    ax.legend()
    plt.show()

    fig_2, ax_2 = plt.subplots(figsize = (10,10))
    ax_2.errorbar(s_full, iq_full, yerr = err_full, fmt= 'o', markersize=3, ecolor="lightgray", label="Experiment")
    ax_2.set_yscale("log")

    ax_2.plot(s_sim, iq_sim, zorder=3, lw=3, label=f_name)
    ax_2.set_ylabel("i(q)")
    ax_2.set_xlabel("s")
    ax_2.set_title("Simulated SAXS fit with Experiment - full")
    ax.legend()
    plt.show()

##########------ MAIN

def main():
    #Find all simulated SAXS curves and place in dataframe
    lent, sangle, intense, file_name, weights_df, wlist = simulated_curves(sim_path)

    #Perform the weighted average SAXS curve
    s_weighted, iq_weighted = average_curve(sangle, intense, weights_df)

    #Create the SAXS curves from experimental data for plotting
    angletrun, intensetrun, errtrun, anglefull, intensefull, errfull = experimental_curve(path_exp_file, lent)

    #Plot experiment vs weighted average simulation
    plot_compare(s_weighted, iq_weighted, angletrun, intensetrun, errtrun, anglefull, intensefull, errfull, file_name)

main()
print("Breakpt")