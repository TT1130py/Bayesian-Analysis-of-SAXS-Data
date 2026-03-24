import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

path_exp_file = "/home/malab/Desktop/BioEn-master/examples/scattering/files/experimental_data"
sim_path = "/home/malab/Documents/sim_curves_test.txt"

def simulated_curves(sim_file):
    sim_pd = pd.DataFrame(pd.read_csv(sim_file, header=None))

    curves = len(sim_pd)

    file = sim_pd.iloc[0,0]
    read_file = pd.DataFrame(pd.read_csv(f"{file}", sep=r"\s+", skiprows=5))
    s_sim = read_file.iloc[:,0] * 10
    iq_sim = read_file.iloc[:,1]

    length = s_sim.iloc[-1]

    return length, s_sim, iq_sim, file

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
#for i in range(curves-1):


len, sangle, intense, file_name = simulated_curves(sim_path)
angletrun, intensetrun, errtrun, anglefull, intensefull, errfull = experimental_curve(path_exp_file, len)
plot_compare(sangle, intense, angletrun, intensetrun, errtrun, anglefull, intensefull, errfull, file_name)
print("Breakpt")