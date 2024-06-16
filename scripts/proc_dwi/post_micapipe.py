"""
Converts the micapipe SC to labeled CSV files
"""

import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


PROJECT_DIR = os.environ.get("PROJECT_DIR")
PNC_PROJECT_DIR = os.environ.get("PNC_PROJECT_DIR")
IMAGEN_PROJECT_DIR = os.environ.get("IMAGEN_PROJECT_DIR")
CODE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def get_output_dir(dataset, ses):
    """
    Determines output directory based on dataset
    """
    if dataset == "pnc":
        OUTPUT_DIR = os.path.join(PNC_PROJECT_DIR, "output")
    elif dataset == "imagen":
        OUTPUT_DIR = os.path.join(IMAGEN_PROJECT_DIR, "output", ses)
    else:
        OUTPUT_DIR = os.path.join(PROJECT_DIR, "output", dataset)
    return OUTPUT_DIR


def label_sc(
    dataset,
    participant_label,
    ses,
    parc,
    measure,
    hemi=None,
    exc_subcortex=True,
    exc_interhemispheric=False,
):
    """
    Labels, normalizes and filters the SC output from micapipe

    Parameters
    ----------
    dataset: (str)
    participant_label: (str)
    ses: (str | None)
    parc: (str)
    measure: (str)
        'strength' or 'length'
    hemi: (str | None)
        limit the SC to hemisphere, if indicated
    exc_subcortex: (bool)
    exc_interhemispheric: (bool)

    Returns
    -------
    sym_SC: (pd.DataFrame)
    """
    OUTPUT_DIR = get_output_dir(dataset, ses)
    # specify subject micapipe directory
    sub_micapipe_dir = os.path.join(OUTPUT_DIR, "micapipe", participant_label)
    if not os.path.exists(sub_micapipe_dir):
        sub_micapipe_dir = os.path.join(OUTPUT_DIR, "micapipe_bak", participant_label)
    # determine SC file prefixes
    if dataset == "micamics":
        SC_file_prefix = os.path.join(
            sub_micapipe_dir,
            "ses-01",
            "dwi",
            f'{participant_label}_ses-01_space-dwinative_atlas-{parc.replace("-", "")}_desc',
        )
    SC_file_prefix = os.path.join(
        sub_micapipe_dir,
        "dwi",
        "connectomes",
        f"{participant_label}_space-dwi_atlas-{parc}_desc-iFOD2-10M-SIFT2_full",
    )
    if measure == "strength":
        if dataset == "micamics":
            SC_file_path = SC_file_prefix + "-sc"
        else:
            SC_file_path = SC_file_prefix + "-connectome"
    elif measure == "length":
        if dataset == "micamics":
            SC_file_path = SC_file_prefix + "-edgeLength"
        else:
            SC_file_path = SC_file_prefix + "-edgeLengths"
    if dataset == "micamics":
        full_SC = np.loadtxt(SC_file_path + ".txt", delimiter=",")
    else:
        full_SC = np.loadtxt(SC_file_path + ".txt")
    # label it and select indicated parcels
    # load parcel info for subcortex and the given parcellation
    parcs_dir = os.path.join(CODE_DIR, "src", "parcellations")
    lut_sctx_mics = pd.read_csv(
        os.path.join(parcs_dir, "lut_subcortical-cerebellum_mics.csv")
    )
    lut_parc_mics = pd.read_csv(os.path.join(parcs_dir, f"lut_{parc}_mics.csv"))
    # order them similar to  micapipe/functions/connectome_slicer.R
    # (which is the order of parcels in full_SC)
    lut_full = pd.concat([lut_parc_mics, lut_sctx_mics], axis=0).sort_values(by="mics")
    # specify parcels to exclude
    ## cerebellum (mics 100-1000)
    if dataset == "micamics":
        # in micamics the subcortical+cortical SC is published and cerebellum parcels are excluded
        # from the SC file => they should also be excluded from lut_full (which later determines
        # labels of full_SC dataframe)
        lut_full = lut_full.loc[
            ~((lut_full["mics"] >= 100) & (lut_full["mics"] < 1000))
        ]
        exc_parcels = []
    else:
        # in other datasets cerebellum parcels exist in full SC but will be excluded in this function
        exc_parcels = lut_full.loc[
            (lut_full["mics"] >= 100) & (lut_full["mics"] < 1000), "mics"
        ].values.tolist()
    ## midline
    exc_parcels += [1000, 2000]
    ## exclude subcortex (mics 10-100) if indicated
    if exc_subcortex:
        exc_parcels += lut_full.loc[(lut_full["mics"]) < 100, "mics"].values.tolist()
    ## exclude the cortical and subcortical parcels from the other hemisphere
    all_hemi_parcs = {
        "L": lut_full.loc[
            ((lut_full["mics"] > 1000) & (lut_full["mics"] < 2000))
            | (lut_full["mics"] < 49),
            "mics",
        ].values.tolist(),
        "R": lut_full.loc[
            (lut_full["mics"] > 2000)
            | ((lut_full["mics"] >= 49) & (lut_full["mics"] < 100)),
            "mics",
        ].values.tolist(),
    }
    if hemi == "L":
        exc_parcels += all_hemi_parcs["R"]
    elif hemi == "R":
        exc_parcels += all_hemi_parcs["L"]
    exc_parcels = list(set(exc_parcels))  # drops duplicates
    # exclude specified parcles and label
    full_SC = pd.DataFrame(full_SC, index=lut_full["mics"], columns=lut_full["mics"])
    SC = full_SC.drop(index=exc_parcels, columns=exc_parcels)
    labels = lut_full.loc[lut_full["mics"].isin(SC.index), "label"].values
    SC.index = SC.columns = labels
    sym_SC = SC + SC.T
    sym_SC.values[np.diag_indices_from(sym_SC)] = 0
    if hemi is None:
        hemi_parcs = {
            "L": sym_SC.index.intersection(
                lut_full.loc[lut_full["mics"].isin(all_hemi_parcs["L"]), "label"]
            ),
            "R": sym_SC.index.intersection(
                lut_full.loc[lut_full["mics"].isin(all_hemi_parcs["R"]), "label"]
            ),
        }
        if not exc_subcortex:
            # order L -> R (instead of ctx -> sctx)
            ordered_parcs = hemi_parcs["L"].tolist() + hemi_parcs["R"].tolist()
            sym_SC = sym_SC.loc[ordered_parcs, ordered_parcs]
        if exc_interhemispheric:
            sym_SC.loc[hemi_parcs["L"], hemi_parcs["R"]] = 0
            sym_SC.loc[hemi_parcs["R"], hemi_parcs["L"]] = 0
    return sym_SC


def normalize_sc_strength(SC, max1=False, mean001=True, threshold_q=1):
    """
    Normalizes raw SC matrix (streamlines count between pairs of
    parcels)

    Parameters
    ---------
    SC: (pd.DataFrame)
    max1: (bool)
        normalize to 0-1
    mean001: (bool)
        normalize mean of SC to 0.01
    threshold_q: (float)
        keep the top threshold_q edges and set the rest to zero

    Returns
    -------
    norm_SC: (pd.DataFrame)
    """
    norm_SC = SC.copy()
    if threshold_q < 1:
        norm_SC.values[
            norm_SC.values
            <= np.quantile(
                norm_SC.values[np.tril_indices_from(norm_SC)], 1 - threshold_q
            )
        ] = 0
    if mean001:
        norm_SC /= norm_SC.values.mean() * 100
    if max1:
        norm_SC /= norm_SC.values.max()
    return norm_SC


def post_micapipe(
    dataset,
    participant_label,
    ses=None,
    parc="schaefer-100",
    hemi=None,
    exc_subcortex=True,
    exc_interhemispheric=False,
    max1=False,
    mean001=True,
    threshold_q=1,
):
    """
    Applies transformations to the micapipe SC output to make it
    ready for simulations

    Parameters
    ----------
    dataset: (str)
    participant_label: (str)
    ses: (str | None)
    parc: (str)
    hemi: (str)
        limit the SC to hemisphere
    exc_subcortex: (bool)
    exc_interhemispheric: (bool)
    max1: (bool)
        normalize to 0-1
    mean001: (bool)
        normalize mean of SC to 0.01
    threshold_q: (float)
        keep the top threshold_q % edges and set the rest to zero
    """
    OUTPUT_DIR = get_output_dir(dataset, ses)
    # specify subject directory
    sub_out_dir = os.path.join(OUTPUT_DIR, "SC", participant_label)
    os.makedirs(sub_out_dir, exist_ok=True)
    # label and reorder SC files and select parcels
    SC_strength = label_sc(
        dataset,
        participant_label,
        ses,
        parc,
        "strength",
        hemi,
        exc_subcortex=exc_subcortex,
        exc_interhemispheric=exc_interhemispheric,
    )
    SC_distance = label_sc(
        dataset,
        participant_label,
        ses,
        parc,
        "length",
        hemi,
        exc_subcortex=exc_subcortex,
        exc_interhemispheric=exc_interhemispheric,
    )
    # normalization of sc trength
    SC_strength = normalize_sc_strength(
        SC_strength,
        max1,
        mean001,
        threshold_q,
    )
    # set distance of absent connections to 0
    if threshold_q < 1:
        SC_distance.values[SC_strength.values == 0] = 0
    # save
    out_prefix = "ctx" if exc_subcortex else "sctx_ctx"
    out_prefix += f"_parc-{parc}"
    if hemi:
        out_prefix += f"_hemi-{hemi}"
    if exc_interhemispheric:
        out_prefix += "_exc-inter"
    if mean001:
        out_prefix += "_mean001"
    if max1:
        out_prefix += "_max1"
    out_prefix += f'_thresh-{str(threshold_q).replace(".", "")}'
    out_prefix = os.path.join(sub_out_dir, out_prefix)
    # save as csv, txt and png
    SC_strength.to_csv(out_prefix + "_desc-strength.csv.gz", compression="gzip")
    SC_distance.to_csv(out_prefix + "_desc-length.csv.gz", compression="gzip")
    np.savetxt(out_prefix + "_desc-strength.txt", SC_strength.values)
    np.savetxt(out_prefix + "_desc-length.txt", SC_distance.values)
    fig, ax = plt.subplots(1, figsize=(8, 7))
    sns.heatmap(SC_strength, cbar=True, ax=ax)
    fig.tight_layout()
    fig.savefig(out_prefix + "_desc-strength.png", dpi=80)
    plt.close(fig)


if __name__ == "__main__":
    dataset = sys.argv[1]
    participant_label = sys.argv[2]
    if dataset == "imagen":
        # we expect sub-xyz_ses-AB as the second argument
        participant_label, ses = participant_label.split("_")
        ses = ses.replace("ses-", "")
    else:
        ses = None
    parcellations = ["schaefer-100", "schaefer-200"]
    for parc in parcellations:
        post_micapipe(
            dataset,
            participant_label,
            ses,
            parc,
            max1=False,
            mean001=True,
        )
