#####----- Compare SAXS curves between multiple sets of reweighed structures

#####----- IMPORTS
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from natsort import natsorted
import seaborn as sns
from scipy.stats import linregress
import matplotlib.lines as mlines

#####----- PATHS, ARGUMENTS, REQUIREMENTS
experimental_data = "/users/t/j/tjaglal/experimental_data/SASDLU4.dat"
structure_path = "/users/t/j/tjaglal/structures"
output_path = ''

## Edit here and/or add another subset for every weight file used
post_weights_1 = ""
grid_sum_1 = "" #this is grid sum
save_path_1 = ""
best_dro_1 = ""
best_r0_1 = ""

wam = ""

## Add new variables to list
dro_list = [best_dro_1]
r0_list = [best_r0_1]
post_weights_list = [post_weights_1]
grid_paths_list = [grid_sum_1]
save_paths_list = [save_path_1]

#####----- FUNCTIONS
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

def rg_guinier(s, iq, fit_points):
    valid_idx = iq > 0
    s_valid = s[valid_idx]
    iq_valid = iq[valid_idx]

    s_low = s_valid[:fit_points]
    iq_low = iq_valid[:fit_points]

    x = s_low**2
    y = np.log(iq_low)

    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    if slope >= 0:
        print("Positive slope")
        return np.nan, r_value**2
    rg = np.sqrt(-3 * slope)

    return rg, r_value**2
def match_files(sim_file, save_path, dro, r0):
    sim_pd = pd.read_csv(sim_file, sep='\\t', header=0)
    sim_pd["PDB_Name"] = sim_pd["PDB_Name"].str.replace('.pdb', '', regex=False)

    grid_df = pd.read_csv(os.path.join(save_path, "grid_full.txt"), sep='\s+', header=None,
                          names=['index', 'dro', 'r0'])
    gp = grid_df.loc[(grid_df['dro'] == dro) & (grid_df['r0'] == r0), 'index'].iloc[0]

    compiled_saxs = np.genfromtxt("{}/compiled_GPs/GP{}_all_saxs.txt".format(save_path, str(gp)))
    compiled_df = pd.DataFrame(compiled_saxs).reset_index(drop=True)
    iq_sim_matrix = pd.DataFrame(compiled_df.drop(columns=[0]).values)

    true_pdb_order = []

    search_pattern = os.path.join(save_path, "mm*", f"GP{gp}", "calc_saxs.txt")
    sorted_files = natsorted(glob.glob(search_pattern))

    for file in sorted_files:
        parts = file.split(os.sep)
        frac_folder = [p for p in parts if p.startswith('mm') and '_' in p][0]

        parent_folder = frac_folder.split('_')[0]
        struct_frac_path = os.path.join(structure_path, parent_folder, frac_folder)
        pdbs = glob.glob(os.path.join(struct_frac_path, "mm*.pdb"))

        pdbs.sort()

        for pdb in pdbs:
            true_pdb_order.append(os.path.basename(pdb).replace('.pdb', ''))

    pdb_names = pd.Series(true_pdb_order)

    s_full = np.genfromtxt("{}/mm016_100/qvals.txt".format(save_path), skip_header=0)
    s_values = s_full[:iq_sim_matrix.shape[1]]

    return s_values, iq_sim_matrix, None, pdb_names

def VACC_average_curve(sim_file, pdb_names, s_val, iq_val):
    sim_pd = pd.read_csv(sim_file, sep='\\t', header=0)
    sim_pd["PDB_Name"] = sim_pd["PDB_Name"].str.replace('.pdb','',regex=False)

    weight_map = dict(zip(sim_pd['PDB_Name'], sim_pd['1']))
    ordered_weights = np.array([weight_map.get(name, 0.0) for name in pdb_names])
    iq_array = np.array(iq_val)

    prior_iq = np.mean(iq_array, axis=1)

    weighted_matrix = iq_array * ordered_weights[:, np.newaxis]
    avg_iq = np.sum(weighted_matrix, axis=0)

    sim_merge = pd.concat([pd.DataFrame(s_val), pd.DataFrame(avg_iq), pd.DataFrame(prior_iq)], axis=1)

    return len(sim_merge), sim_merge.iloc[:, 0], sim_merge.iloc[:, 1], sim_merge.iloc[:, 2]
    print("Breakpt")

def plot_curves(s_sim_arr, iq_sim_arr, iq_prior_arr, s, iq, err, post_rg_arr, pri_rg_arr, exp_rg):
    colors_post = sns.color_palette("flare", n_colors=wam)
    colors_pri = sns.color_palette("crest", n_colors=wam)

    fig, ax = plt.subplots(figsize = (10,10))
    ax.errorbar(s, iq, yerr = err, fmt= 'o', markersize=3, ecolor="lightgray", label="Experiment")
    ax.set_yscale("log")

    for i in range(wam):
        scale = np.sum(iq * iq_sim_arr[i]) / np.sum(iq_sim_arr[i] ** 2)
        scaled_iq_sim = iq_sim_arr[i] * scale

        ax.plot(s_sim_arr[i], scaled_iq_sim, lw=3, label=fr"$\delta \rho$: {dro_list[i]} | $r_0$: {r0_list[i]}", color=colors_post[i])
    ax.set_ylabel("i(q)")
    ax.set_xlabel("s")
    ax.set_title("Posterior SAXS curves with experiment")

    leg_post = ax.legend(loc="upper right")
    ax.add_artist(leg_post)

    rg_handles_post = [mlines.Line2D([], [], color='none', label=f"Exp Rg: {exp_rg:.2f} Å")]
    for i in range(wam):
        label_text = fr"$\delta \rho$: {dro_list[i]} | $r_0$: {r0_list[i]} Rg: {post_rg_arr[i]:.2f} Å"
        rg_handles_post.append(mlines.Line2D([], [], color='none', label=label_text))

    ax.legend(handles=rg_handles_post, loc="lower left", title="Radius of Gyration ($R_g$)", handlelength=0,
              handletextpad=0)

    fig_2, ax_2 = plt.subplots(figsize=(10, 10))
    ax_2.errorbar(s, iq, yerr=err, fmt='o', markersize=3, ecolor="lightgray", label="Experiment")
    ax_2.set_yscale("log")

    for i in range(wam):
        scale_p = np.sum(iq * iq_prior_arr[i]) / np.sum(iq_prior_arr[i] ** 2)
        scaled_iq_prior = iq_prior_arr[i] * scale_p

        ax_2.plot(s_sim_arr[i], scaled_iq_prior, lw=3, label=fr"$\delta \rho$: {dro_list[i]} | $r_0$: {r0_list[i]}", color=colors_pri[i])
    ax_2.set_ylabel("i(q)")
    ax_2.set_xlabel("s")
    ax_2.set_title("Prior SAXS curves with experiment")

    leg_pri = ax_2.legend(loc="upper right")
    ax_2.add_artist(leg_pri)

    # 2nd Legend (Rg Values)
    rg_handles_pri = [mlines.Line2D([], [], color='none', label=f"Exp Rg: {exp_rg:.2f} Å")]
    for i in range(wam):
        label_text = fr"$\delta \rho$: {dro_list[i]} | $r_0$: {r0_list[i]} Rg: {pri_rg_arr[i]:.2f} Å"
        rg_handles_pri.append(mlines.Line2D([], [], color='none', label=label_text))

    ax_2.legend(handles=rg_handles_pri, loc="lower left", title="Radius of Gyration ($R_g$)", handlelength=0,
                handletextpad=0)

    return fig, fig_2

def main():
    s_arr = []
    wiq_arr = []
    piq_arr = []
    pos_rg_arr = []
    pri_rg_arr = []

    lent = 1000
    #Will be overwritten for real plotting
    angletrun, intensetrun, errtrun = experimental_curve(experimental_data, lent)
    exp_rg, exp_r2 = rg_guinier(angletrun, intensetrun, 15)

    for i in range(wam):
        s, concat_merge, weights, f_name = match_files(post_weights_list[i], save_paths_list[i], dro_list[i], r0_list[i])
        lent, s_weighted, iq_weighted, iq_prior = VACC_average_curve(post_weights_list[i], f_name, s, concat_merge)
        s_arr.append(s_weighted)
        wiq_arr.append(iq_weighted)
        piq_arr.append(iq_prior)

        post_rg, post_r2 = rg_guinier(s_weighted, iq_weighted, 15)
        prior_rg, prior_r2 = rg_guinier(s_weighted, iq_prior, 15)
        pos_rg_arr.append(post_rg)
        pri_rg_arr.append(prior_rg)

    angletrun, intensetrun, errtrun = experimental_curve(experimental_data, lent)
    post_figure, pri_figure = plot_curves(s_arr, wiq_arr, piq_arr, angletrun, intensetrun, errtrun, exp_rg, pos_rg_arr, pri_rg_arr)

    post_figure.savefig(f"{output_path}/posterior_saxs_curves.png", dpi=300, bbox_inches='tight')
    pri_figure.savefig(f"{output_path}/prior_saxs_curves.png", dpi=300, bbox_inches='tight')

if __name__ == '__main__':
    main()