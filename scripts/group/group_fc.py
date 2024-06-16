"""
Calculates group-representative FC matrix of the given dataset and subjects
"""

import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

PROJECT_DIR = os.environ.get("PROJECT_DIR")
PNC_PROJECT_DIR = os.environ.get("PNC_PROJECT_DIR")
CODE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(CODE_DIR)
from proc_rest.post_fmriprep import *


def calculate_group_fc(
    dataset, fc_file_prefix, group_name, subjects=[], downsample_fcd=True
):
    """
    Calculates group-averaged FC

    Parameters
    ----------
    dataset: (str)
    fc_file_prefix: (str)
        e.g. ctx_parc-schaefer-100_hemi-LR_highpass-013_lowpass-none_desc-FC
    subjects: (list)
        if empty gets all the subjects from the dataset inside FC directory
    downsample_fcd: (bool)
        if true downsamples FCD by np.sort(concat_fcd)[::len(subjects)]
    """
    if dataset == "pnc":
        OUTPUT_DIR = os.path.join(PNC_PROJECT_DIR, "output")
    else:
        raise NotImplementedError("Only PNC is supported")
    fc_dir = os.path.join(OUTPUT_DIR, "FC")
    # get all the subjects if an empty list is provided
    if group_name == "all":
        for subdir in os.listdir(fc_dir):
            if "sub-" in subdir:
                subjects.append(subdir)
        subjects = sorted(subjects)
    else:
        if len(subjects) < 2:
            raise ValueError("Group must contain at least 2 subjects")
    # read subject FCs and z-transform them
    z_FCs = []
    included_subjects = []
    FCD_trils = []
    for subject in subjects:
        sub_fc_file = os.path.join(fc_dir, subject, fc_file_prefix + ".csv.gz")
        sub_fcd_file = os.path.join(fc_dir, subject, fc_file_prefix + "Dtril.txt")
        if os.path.exists(sub_fc_file):
            FC = pd.read_csv(sub_fc_file, compression="gzip", index_col=0)
            parcels = FC.index
            FC = FC.values
            # zero out correlations of 1 (to avoid division by 0)
            FC[np.isclose(FC, 1)] = 0
            # Fisher's z-transformation
            FC = np.arctanh(FC)
            # # zero out NaNs and inf
            # FC[np.isnan(FC) | np.isinf(FC)] = 0
            z_FCs.append(FC[:, :, np.newaxis])
            included_subjects.append(subject)
            ## FCD
            FCD_trils.append(np.loadtxt(sub_fcd_file))
    # concatenate and average across subjects
    z_FCs = np.concatenate(z_FCs, axis=2)
    group_z_FC = z_FCs.mean(axis=2)
    # inverse z to r and add back diagonal of 1s
    group_r_FC = np.tanh(group_z_FC)
    group_r_FC[np.diag_indices_from(group_r_FC)] = 1.0
    # convert to dataframe and save as csv, png and txt (for lower triangle)
    group_r_FC = pd.DataFrame(group_r_FC, index=parcels, columns=parcels)
    out_dir = os.path.join(fc_dir, f"group-{group_name}")
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "subjects.txt"), "w") as f:
        # store the list of subjects (important in case one subject misses
        # a particular session but not the others)
        f.write("\n".join(included_subjects))
    out_prefix = os.path.join(out_dir, fc_file_prefix)
    group_r_FC.to_csv(out_prefix + ".csv.gz", compression="gzip")
    fig, ax = plt.subplots(1, figsize=(4, 4))
    sns.heatmap(group_r_FC.values, cmap="RdBu_r", vmin=-1, vmax=1, cbar=False, ax=ax)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_prefix + ".png", dpi=80)
    plt.close(fig)
    fc_tril = group_r_FC.values[np.tril_indices_from(group_r_FC, -1)]
    fc_tril = fc_tril[~np.isnan(fc_tril)]
    np.savetxt(out_prefix + "tril.txt", fc_tril, fmt="%.8f")
    ## FCD
    pooled_FCD_tril = np.concatenate(FCD_trils)
    if downsample_fcd:
        pooled_FCD_tril = np.sort(pooled_FCD_tril)[:: len(subjects)]
    np.savetxt(out_prefix + "Dtril.txt", pooled_FCD_tril, fmt="%.8f")


if __name__ == "__main__":
    dataset = sys.argv[1]
    fc_file_prefix = sys.argv[2]
    if len(sys.argv) > 3:
        group_name = sys.argv[3]
        subjects = sys.argv[4:]
    else:
        group_name = "all"
        subjects = []
    subjects = sys.argv[3:]  # if no subjects are provided the list will be empty
    calculate_group_fc(dataset, fc_file_prefix, group_name, subjects)
