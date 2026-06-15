import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns

#####----- PATHS, ARGUMENTS, REQUIREMENTS
path_exp_file = ""
search_path = ""
save_path = ""

files = []
#####----- FUNCTIONs

def experimental_curve(path_exp_file, sim_length, save_path):
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
    #save_path_1 = "{}/experiment.png".format(save_path)
    #plt.savefig(save_path, dpi=300)

    return s_trun, iq_trun, err_trun

def find_curves(files, search_path):
    q = []
    iq = []
    for i in range(len(files)):
        pdb = pd.read_csv(f"{search_path}/{files[i]}.out", sep=r"\s+", skiprows=5, header=0)
        pdb.columns = pdb.columns.str.lstrip("#")

        q_np = pdb["q"].to_numpy()
        iq_np = pdb["I"].to_numpy()

        q.append(q_np)
        iq.append(iq_np)

    lent = len(q_np[0])
    return lent, q, iq

def plot_curves(s, iq, err, files, q_arr, iq_arr):
    colors = sns.color_palette("crest", n_colors=len(files))

    fig, ax = plt.subplots(figsize = (10,10))
    ax.errorbar(s, iq, yerr = err, fmt= 'o', markersize=3, ecolor="lightgray", label="Experiment")
    ax.set_yscale("log")

    for i in range(len(files)):
        ax.plot(q_arr[0], iq_arr[i], zorder= i+1, lw=3, label= files[i], color=colors[i])
    ax.set_ylabel("i(q)")
    ax.set_xlabel("s")
    ax.set_title("SAXS curves vs Experiment")
    ax.legend()

    save_path_2 = "{}/ind_curves.png".format(save_path)
    fig.savefig(save_path_2, dpi=300)

def main():
    length, q_files, iq_files = find_curves(files, search_path)
    s_exp, iq_exp, err_exp = experimental_curve(path_exp_file, length, save_path)
    plot_curves(s_exp, iq_exp, err_exp, files, q_files, iq_files)


