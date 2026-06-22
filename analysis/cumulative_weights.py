import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse
import yaml

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

    return fig

#####---------- MAIN
#####----- Initialize CLI arguments and configuration
parser = argparse.ArgumentParser()

parser.add_argument("--config", type=str, default="config.yaml", help="Path to main YAML file")
args = parser.parse_parser_args() if hasattr(parser, 'parse_parser_args') else parser.parse_args()

with open(args.config, "r") as f:
    master_config = yaml.safe_load(f)

cw_config = master_config.get("cumulative_weights", {})


reweighted_path = cw_config.get["reweighted_path", ""]
save_path = cw_config.get["save_path", ""]

pos_w, pri_w, index = cumulative_weights_all(reweighted_path)
fig = plot_cumulative(pos_w, pri_w, index)
fig.savefig(os.path.join(save_path, "cumul_weights_r0r_t2000.png"))
