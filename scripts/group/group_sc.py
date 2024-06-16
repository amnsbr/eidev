import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


PROJECT_DIR = os.environ.get("PROJECT_DIR")
PNC_PROJECT_DIR = os.environ.get("PNC_PROJECT_DIR")
CODE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(CODE_DIR)
from proc_dwi import post_micapipe


def calculate_group_sc(
    dataset,
    parc,
    hemi,
    group_name,
    subjects=[],
    ses=None,
    approach="median",
    threshold=1.0,
    max1=False,
    mean001=True,
    exc_subcortex=True,
    exc_interhemispheric=False,
    overwrite=False,
):
    """
    Calculates group-averaged SC

    Parameters
    ----------
    dataset: (str)
    parc: (str)
    hemi: (str | None)
    subjects: (list)
        if empty gets all the subjects from the dataset inside SC directory
    ses: (str | None)
    approach: (str)
        - 'mean'
        - 'median'
    threshold: (float)
    max1: (bool)
        normalize by division to max
    mean001: (bool)
        normalize by division to (mean * 100)
    exc_subcortex: (bool)
    exc_interhemispheric: (bool)
    overwrite: (bool)
    """
    # specify output paths
    OUTPUT_DIR = post_micapipe.get_output_dir(dataset, ses)
    out_dir = os.path.join(OUTPUT_DIR, "SC", f"group-{group_name}")
    out_prefix = "ctx" if exc_subcortex else "sctx_ctx"
    out_prefix += f"_parc-{parc}"
    out_prefix += f"_approach-{approach}"
    if hemi:
        out_prefix += f"_hemi-{hemi}"
    if exc_interhemispheric:
        out_prefix += "_exc-inter"
    if mean001:
        out_prefix += "_mean001"
    if max1:
        out_prefix += "_max1"
    out_prefix = os.path.join(out_dir, out_prefix)
    # do not recreate it if it exists
    if not overwrite:
        sc_strength_path = out_prefix + "_desc-strength.csv.gz"
        sc_length_path = out_prefix + "_desc-length.csv.gz"
        if os.path.exists(sc_strength_path) & os.path.exists(sc_length_path):
            group_sc_strength = pd.read_csv(sc_strength_path, index_col=0)
            group_sc_length = pd.read_csv(sc_length_path, index_col=0)
            return group_sc_strength, group_sc_length

    # get all the subjects if an empty list is provided
    SC_DIR = os.path.join(OUTPUT_DIR, "SC")
    if group_name == "all":
        subjects = []
        for subdir in os.listdir(SC_DIR):
            if not ("group-" in subdir):
                subjects.append(subdir)
        subjects = sorted(subjects)
    else:
        if len(subjects) < 2:
            raise ValueError("Group must contain at least 2 subjects")
    # load SC strength and length for each subject (as np.ndarray)
    sc_strengths = []
    sc_lengths = []
    included_subjects = []
    for subject in subjects:
        try:
            sc_strength_df = post_micapipe.label_sc(
                dataset,
                subject,
                ses,
                parc,
                "strength",
                hemi,
                exc_subcortex=exc_subcortex,
                exc_interhemispheric=exc_interhemispheric,
            )
        except FileNotFoundError:
            print(f"Connectome not found for {subject}")
            continue
        sc_strengths.append(sc_strength_df.values[:, :, np.newaxis])
        sc_length_df = post_micapipe.label_sc(
            dataset,
            subject,
            ses,
            parc,
            "length",
            hemi,
            exc_subcortex=exc_subcortex,
            exc_interhemispheric=exc_interhemispheric,
        )
        sc_lengths.append(sc_length_df.values[:, :, np.newaxis])
        included_subjects.append(subject)
    parcels = sc_strength_df.index
    # aggregate SC strength and length across subjects
    sc_strengths = np.concatenate(sc_strengths, axis=2)
    sc_lengths = np.concatenate(sc_lengths, axis=2)
    if approach == "mean":
        group_sc_strength = np.nanmean(sc_strengths, axis=2)
        group_sc_length = np.nanmean(sc_lengths, axis=2)
    elif approach == "median":
        group_sc_strength = np.nanmedian(sc_strengths, axis=2)
        group_sc_length = np.nanmedian(sc_lengths, axis=2)
    # convert to labeled df
    group_sc_strength = pd.DataFrame(group_sc_strength, index=parcels, columns=parcels)
    group_sc_length = pd.DataFrame(group_sc_length, index=parcels, columns=parcels)
    # normalize SC strength if indicated
    group_sc_strength = post_micapipe.normalize_sc_strength(
        group_sc_strength, max1=max1, mean001=mean001, threshold_q=threshold
    )
    # set distance of absent connections to 0
    group_sc_length.values[group_sc_strength.values == 0] = 0
    # save
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "subjects.txt"), "w") as f:
        # store the list of subjects (important in case one subject misses
        # a particular session but not the others)
        f.write("\n".join(included_subjects))
    # save as csv, txt and png
    group_sc_strength.to_csv(out_prefix + "_desc-strength.csv.gz", compression="gzip")
    group_sc_length.to_csv(out_prefix + "_desc-length.csv.gz", compression="gzip")
    np.savetxt(out_prefix + "_desc-strength.txt", group_sc_strength.values)
    np.savetxt(out_prefix + "_desc-length.txt", group_sc_length.values)
    fig, ax = plt.subplots(1, figsize=(8, 7))
    sns.heatmap(group_sc_strength, cbar=True, ax=ax)
    fig.tight_layout()
    fig.savefig(out_prefix + "_desc-strength.png", dpi=80)
    plt.close(fig)
    return group_sc_strength, group_sc_length


if __name__ == "__main__":
    dataset = sys.argv[1]
    if len(sys.argv) > 2:
        group_name = sys.argv[2]
        subjects = sys.argv[3:]
    else:
        group_name = "all"
        subjects = []
    parcellations = ["schaefer-100", "schaefer-200"]
    for parc in parcellations:
        calculate_group_sc(dataset, parc, None, group_name, subjects)
