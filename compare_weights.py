import pandas as pd

weights_path_1 = "/home/malab/Downloads/structure_weights_sorted_2025-11-22_0.csv"
path_1_type = "SAXS"
weights_path_2 = "/home/malab/Downloads/Sorted_Results_NMR_Full.csv"
path_2_type = "NMR"
save_path = "/home/malab/Desktop/compare_weights_output" #optional save path

## Functions -------------
def match_weights(w1,w2):
    w1_df = pd.read_csv(weights_path_1, header=None)
    w2_df = pd.read_csv(weights_path_2, header=None)

    #Combining weights for each file
    combined_weights = []
    for i in range(len(w1_df)):
        pdb_file = w1_df.loc[i,0]
        w1_weight = w1_df.loc[i,1]

        mask = w2_df[0] == pdb_file
        w2_rows = w2_df.loc[mask, 1]

        if not w2_rows.empty:
            w2_weight = w2_rows.iloc[0]   #take first match
        else:
            w2_weight = None

        #add weights to new dataframe
        combined_weights.append([pdb_file, w1_weight, w2_weight])

    return combined_weights

def rank_weights(weight_table):
    cw_df = pd.DataFrame(weight_table)

    #take average weight between two datasets
    averages = []
    for i in range(len(cw_df)):
        avg = (cw_df.loc[i,1] + cw_df.loc[i,2]) / 2
        averages.append(avg)
    avg_df = pd.DataFrame(averages)
    cw_df['Average weight'] = avg_df

    #reorganize by highest averages
    ranked_df = cw_df.sort_values(by='Average weight', ascending=False)
    name_ranked_df = ranked_df.rename(columns={0: "PDB File", 1: path_1_type, 2: path_2_type})

    return name_ranked_df

## Main ----------------
df = match_weights(weights_path_1, weights_path_2)
final = rank_weights(df)
print("Breakpt")

##Optional save file as csv
final.to_csv("{}/compared_weights.csv".format(save_path))