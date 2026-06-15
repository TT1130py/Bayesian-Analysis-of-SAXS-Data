import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


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
    save_path_1 = "{}/experiment.png".format(save_path)
    plt.savefig(save_path_1, dpi=300)

    return s_trun, iq_trun, err_trun