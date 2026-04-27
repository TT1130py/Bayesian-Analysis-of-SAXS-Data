import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

######--------- PATHS

reweighted_path = "/home/malab/Desktop/compare_weights_info/structure_weights_sorted_2026-04-24.txt"
save_path = "/home/malab/Desktop/compare_weights_info" #optional

#####--------- FUNCTIONS

def cumulative_weights_all(pos_path):
    posterior = pd.read_csv(pos_path, sep='\t')

    #Extract just weights
    w_posterior = posterior.iloc[:,1]

    #Define prior weights
    frames = len(w_posterior)
    w_prior = np.ones(frames) / frames
    w_index = np.arange(frames)

    #Desceding sort for plotting purposes
    sorted_w_posterior = np.sort(w_posterior)[::-1]
    sorted_w_prior = np.sort(w_prior)[::-1]

    #Cumulative sums
    posterior_cumulative = np.cumsum(sorted_w_posterior)
    prior_cumulative = np.cumsum(sorted_w_prior)

    return posterior_cumulative, prior_cumulative, w_index

def plot_cumulative(pos, pri, idx):
    fig, ax = plt.subplots(figsize=(10,10))

    ax.plot(idx, pri, label="Prior", color='red')
    ax.plot(idx, pos, label="Posterior", color='blue')

    ax.set_xlabel("Weight index")
    ax.set_ylabel("Cumulative weight")
    ax.set_title("Cumulative weights of ordered w for posterior and prior")
    ax.legend()
    plt.tight_layout()
    plt.show()

    fig.savefig(os.path.join(save_path, "cumul_weights.png"))

#####---------- MAIN
pos_w, pri_w, index = cumulative_weights_all(reweighted_path)
plot_cumulative(pos_w, pri_w, index)