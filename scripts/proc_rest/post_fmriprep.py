"""
Performs additional cleaning and parcellation on the
output of fmriprep
"""

import os, sys
import json
import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt
import seaborn as sns
import nilearn.image, nilearn.surface, nilearn.signal
from nilearn.input_data import NiftiLabelsMasker

PROJECT_DIR = os.environ.get("PROJECT_DIR")
PNC_PROJECT_DIR = os.environ.get("PNC_PROJECT_DIR")
IMAGEN_PROJECT_DIR = os.environ.get("IMAGEN_PROJECT_DIR")
CODE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
SRC_DIR = os.path.join(CODE_DIR, "src")
sys.path.append(CODE_DIR)
from utils import transform

# windows size ~ 30s
WINDOW_SIZE = {
    "pnc": 10,  # TR = 3
    "imagen": 14,  # TR = 2.2
}

# step ~ 5s
STEP = {
    "pnc": 2,
    "imagen": 2,
}

MIN_GOOD_SCAN_DURATION = {  # minutes
    "pnc": 4,
    "imagen": 4,
}

TRs = {  # seconds
    "pnc": 3,
    "imagen": 2.2,
}


def calculate_fcd(
    bold,
    window_size=20,
    step=4,
    drop_edges=True,
    exc_interhemispheric=True,
    hemi_parcs={},
):
    """
    Calculates FC dynamics matrix

    Parameters
    ---------
    bold: (pd.DataFrame) n_regions x n_vols
    window_size: (int)
        must be even, the actual window_size is +1 (including center)
    step: (int)
    drop_edges: (bool)
    exc_interhemispheric: (bool)
    hemi_parcs: (dict)
        parcels of 'L' and 'R' hemisphere
        (only needed when exc_interhemispheric is True)

    Returns
    -------
    fcd_matrix: (np.corrcoef) n_windows x n_windows
    window_fcs: (np.corrcoef) n_regions x n_regions x n_windows
    """
    assert window_size % 2 == 0, "Window size must be even"
    n_regions = bold.shape[0]
    n_vols = bold.shape[1]
    window_fc_trils = []
    window_fcs = []
    # calculate the FC of each window
    if drop_edges:
        first_center = int(window_size / 2)
        last_center = n_vols - 1 - int(window_size / 2)
    else:
        first_center = 0
        last_center = n_vols - 1
    window_center = first_center
    n_windows = 0
    while window_center <= last_center:
        window_start = window_center - int(window_size / 2)
        if window_start < 0:
            window_start = 0
        window_end = window_center + int(window_size / 2)
        if window_end >= n_vols:
            window_end = n_vols - 1
        window_bold = bold.values[:, window_start : window_end + 1]
        window_center += step
        # discard the window if more than 50% of its
        # volumes are motion outliers
        if (window_bold.sum(axis=0) == 0).mean() > 0.5:
            continue
        window_fc = pd.DataFrame(
            np.corrcoef(window_bold), index=bold.index, columns=bold.index
        )
        if exc_interhemispheric:
            window_fc.loc[hemi_parcs["L"], hemi_parcs["R"]] = np.NaN
            window_fc.loc[hemi_parcs["R"], hemi_parcs["L"]] = np.NaN
        window_fcs.append(window_fc.to_numpy()[:, :, np.newaxis])
        window_fc_tril = window_fcs[-1][np.tril_indices(n_regions, -1)]
        window_fc_tril = window_fc_tril[~np.isnan(window_fc_tril), np.newaxis]
        window_fc_trils.append(window_fc_tril)
        n_windows += 1
    window_fcs = np.concatenate(window_fcs, axis=2)
    window_fc_trils = np.concatenate(window_fc_trils, axis=1)
    fcd_matrix = np.corrcoef(window_fc_trils.T)
    return fcd_matrix, window_fcs


def post_fmriprep(
    dataset,
    participant_label,
    ses=None,
    parcellation_name="schaefer-100",
    vols_remove=3,
    high_pass=0.013,
    low_pass=None,
    rms_threshold=0.25,
    fd_threshold=None,
    dvars_threshold=None,
    spike_window=False,
    skip_scrubbing=False,
    min_good_scan_duration=4,
    max_mean_motion={"rms": 0.2},
    save_dynamicFC=False,
):
    """
    Performs additional preprocessing steps after fmriprep, including parcellation,
    detrending, high-pass filtering, confound removal, and motion censoring. Saves
    the parcellated BOLD activity in the postfmriprep directories.

    Parameters
    ---------
    dataset: (str)
    participant_label: (str)
    ses: (str | None)
    parcellation_name: (str)
    vols_remove: (int)
    high_pass: (int | None)
    low_pass: (int | None)
    rms_threshold: (float | None)
    fd_threshold: (float | None)
    dvars_threshold: (float | None)
    spike_window: (bool)
        for motion spikes also remove one volume before and two volumes after
    skin_scrubbing: (bool)
    min_good_scan_duration: (int)
        in minutes; the cutoff to exclude subjects
    max_mean_motion: (dict)
        for each motion parameter the cut-off to exclude subjects
    save_dynamicFC: (bool) not recommended (large files)
    """
    if dataset == "pnc":
        OUTPUT_DIR = os.path.join(PNC_PROJECT_DIR, "output")
    elif dataset == "imagen":
        OUTPUT_DIR = os.path.join(IMAGEN_PROJECT_DIR, "output", ses)
    sub_fmriprep_dir = os.path.join(OUTPUT_DIR, "fmriprep", participant_label)
    sub_func_prefix = os.path.join(
        sub_fmriprep_dir, "func", f"{participant_label}_task-rest"
    )
    sub_out_dir = os.path.join(OUTPUT_DIR, "postfmriprep", participant_label)
    participants_tsv_path = os.path.join(OUTPUT_DIR, "postfmriprep", "participants")
    config_str = (
        f'_highpass-{str(high_pass).replace("0.","").lower()}'
        + f'_lowpass-{str(low_pass).replace("0.","").lower()}'
    )
    if skip_scrubbing:
        config_str += "_noscrub"
        fd_threshold = dvars_threshold = None
    os.makedirs(sub_out_dir, exist_ok=True)
    # Check if cleaned bold already exists
    parc_bold_clean = {}
    out_paths = {}
    all_exist = True
    for structure in ["ctx", "sctx"]:
        parc_bold_clean[structure] = {}
        out_paths[structure] = {}
        for hemi in ["L", "R"]:
            out_path = os.path.join(sub_out_dir, structure)
            if structure == "ctx":
                out_path += f"_parc-{parcellation_name}"
            out_path += f"_hemi-{hemi}" + config_str + "_desc-cleanbold.csv.gz"
            out_paths[structure][hemi] = out_path
            all_exist = all_exist and os.path.exists(out_path)
            if all_exist:
                parc_bold_clean[structure][hemi] = pd.read_csv(out_path, index_col=0)
    if all_exist:
        print("Parcellated clean BOLD already calculated")
        # FC and FCD calculation
        sub_fc_dir = sub_out_dir.replace("postfmriprep", "FC")
        _bold_to_fc_fcd(
            parc_bold_clean,
            sub_fc_dir,
            dataset,
            parcellation_name,
            config_str,
            save_dynamicFC,
        )
        return
    # Parcellate bold signal
    print("parcellating BOLD signal")
    parc_bold = {}
    ## cortex
    bold_fsa = {}
    for hemi in ["L", "R"]:
        bold_fsa[hemi] = nilearn.surface.load_surf_data(
            f"{sub_func_prefix}_hemi-{hemi}_space-fsaverage5_bold.func.gii"
        )
    parc_bold["ctx"] = transform.parcellate_surf(bold_fsa, parcellation_name)
    ## subcortex
    bold_mni = nilearn.image.load_img(
        f"{sub_func_prefix}_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz"
    )
    aseg_dseg_path = f"{sub_func_prefix}_space-MNI152NLin6Asym_desc-aseg_dseg.nii.gz"
    masker = NiftiLabelsMasker(aseg_dseg_path)
    bold_aseg_dseg = pd.DataFrame(
        masker.fit_transform(bold_mni), columns=masker.labels_
    )
    lut_subcortical = pd.read_csv(
        os.path.join(SRC_DIR, "parcellations", "lut_subcortical-cerebellum_mics.csv")
    )
    parc_bold["sctx"] = {}
    for hemi in ["left", "right"]:
        hemi_short = hemi[0].upper()
        subcortical_parcels = lut_subcortical[
            (lut_subcortical["mics"] < 100) & (lut_subcortical["side"] == hemi)
        ].set_index("mics")["label"]
        parc_bold["sctx"][hemi_short] = bold_aseg_dseg.loc[
            :, subcortical_parcels.index
        ].T
        parc_bold["sctx"][hemi_short].index = subcortical_parcels.values
    # get TR from sidecar json
    bold_json_path = f"{sub_func_prefix}_space-MNI152NLin6Asym_desc-preproc_bold.json"
    with open(bold_json_path, "r") as jsonf:
        bold_json = json.load(jsonf)
    t_r = bold_json.get("RepetitionTime", TRs[dataset])
    # select confounds
    confounds_all = pd.read_csv(
        f"{sub_func_prefix}_desc-confounds_timeseries.tsv", sep="\t"
    )
    confound_cols = ["trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z"]
    confound_cols += [c + "_derivative1" for c in confound_cols]  # derivative
    confound_cols += [c + "_power2" for c in confound_cols]  # power2
    confound_cols += ["white_matter", "csf"]
    confounds = confounds_all.iloc[vols_remove:, :].loc[:, confound_cols]
    n_vols = confounds.shape[0]
    # find motion outliers and one volume before and two volumes after them
    std_dvars = confounds_all.iloc[vols_remove:, :].loc[:, "std_dvars"].values
    fd = confounds_all.iloc[vols_remove:, :].loc[:, "framewise_displacement"].values
    rms = confounds_all.iloc[vols_remove:, :].loc[:, "rmsd"].values
    if fd_threshold is None:
        fd_threshold = fd.max()
    if dvars_threshold is None:
        dvars_threshold = std_dvars.max()
    if rms_threshold is None:
        rms_threshold = rms.max()
    motion_outliers = np.where(
        (rms > rms_threshold) | (std_dvars > dvars_threshold) | (fd > fd_threshold)
    )[0]
    if spike_window:
        motion_outliers_ext = np.concatenate(
            [
                motion_outliers,
                motion_outliers - 1,
                motion_outliers + 1,
                motion_outliers + 2,
            ]
        )
        motion_outliers_ext[motion_outliers_ext < 0] = 0
        motion_outliers_ext[motion_outliers_ext > fd.shape[0] - 1] = fd.shape[0] - 1
        motion_outliers_ext = np.unique(motion_outliers_ext)
    else:
        motion_outliers_ext = motion_outliers
    good_vols = list(set(range(n_vols)) - set(motion_outliers_ext))
    ## save motion stats to participants.tsv
    try:
        participants_tsv = pd.read_csv(
            participants_tsv_path + ".tsv", sep="\t", index_col=0
        )
    except:
        participants_tsv = pd.DataFrame(
            columns=[
                "n_motion_outliers",
                "pct_motion_outliers",
                "good_scan_duration",
                "mean_rms",
                "mean_fd",
                "mean_std_dvars",
            ]
        )
    participants_tsv.loc[participant_label, "n_motion_outliers"] = len(
        motion_outliers_ext
    )
    participants_tsv.loc[participant_label, "pct_motion_outliers"] = (
        len(motion_outliers_ext) / n_vols
    )
    participants_tsv.loc[participant_label, "good_scan_duration"] = (
        len(good_vols) * t_r / 60
    )  # approximate
    participants_tsv.loc[participant_label, "mean_rms"] = rms.mean()
    participants_tsv.loc[participant_label, "mean_fd"] = fd.mean()
    participants_tsv.loc[participant_label, "mean_std_dvars"] = std_dvars.mean()
    participants_tsv.loc[participant_label, "n_rms_gte_025"] = (rms > 0.25).sum()
    participants_tsv.to_csv(
        participants_tsv_path + ".tsv", sep="\t", index_label="participant_label"
    )
    if (
        participants_tsv.loc[participant_label, "good_scan_duration"]
        < min_good_scan_duration
    ):
        os.makedirs(os.path.join(OUTPUT_DIR, "excluded"), exist_ok=True)
        with open(
            os.path.join(OUTPUT_DIR, "excluded", participant_label + ".txt"), "a"
        ) as f:
            f.write(
                f"Length of usable scan is less than {min_good_scan_duration} minutes."
            )
    for motion_measure in max_mean_motion:
        if (
            participants_tsv.loc[participant_label, f"mean_{motion_measure}"]
            > max_mean_motion[motion_measure]
        ):
            os.makedirs(os.path.join(OUTPUT_DIR, "excluded"), exist_ok=True)
            with open(
                os.path.join(OUTPUT_DIR, "excluded", participant_label + ".txt"), "a"
            ) as f:
                f.write(
                    f"Mean {motion_measure} greater than {max_mean_motion[motion_measure]}."
                )
    # bold cleaning
    print("cleaning BOLD")
    for structure in ["ctx", "sctx"]:
        for hemi in ["L", "R"]:
            # Removing first volumes, high/low-pass filtering and confound removal
            parc_bold_clean[structure][hemi] = nilearn.signal.clean(
                parc_bold[structure][hemi].values[:, vols_remove:].T,  # remove the first volumes
                detrend=False,
                standardize=False,
                confounds=confounds,
                high_pass=high_pass,
                low_pass=low_pass,
                t_r=t_r,
            ).T
            # relabel and drop NA (midline parcel(s))
            parc_bold_clean[structure][hemi] = pd.DataFrame(
                parc_bold_clean[structure][hemi], index=parc_bold[structure][hemi].index
            ).dropna()
            # scrub motion outliers by z-scoring the no-outlier volumes
            # and setting the signal at outliers to zero
            # (for FC calculation this is equivalent to simply
            # excluding motion outliers, but this approach is preferred
            # as it preserves the temporal structure of the signal - important
            # for FCD calculations)
            parc_bold_clean[structure][hemi].iloc[:, good_vols] = scipy.stats.zscore(
                parc_bold_clean[structure][hemi].iloc[:, good_vols], axis=1
            )
            parc_bold_clean[structure][hemi].iloc[:, motion_outliers_ext] = 0
            parc_bold_clean[structure][hemi].to_csv(
                out_paths[structure][hemi], compression="gzip"
            )
            fig, ax = plt.subplots(1, figsize=(20, 8))
            vmin = min(
                parc_bold_clean[structure][hemi].values.min(),
                -parc_bold_clean["ctx"]["L"].values.max(),
            )
            sns.heatmap(
                parc_bold_clean[structure][hemi].values,
                cmap="RdBu",
                vmin=vmin,
                vmax=-vmin,
                cbar=True,
                ax=ax,
            )
            for spike in motion_outliers_ext:
                ax.axvline(spike + 0.5, linewidth=3, color="black")
            ax.axis("off")
            fig.tight_layout()
            fig.savefig(out_paths[structure][hemi].replace(".csv.gz", ".png"), dpi=80)
            plt.close(fig)
    # FC and FCD calculation
    sub_fc_dir = sub_out_dir.replace("postfmriprep", "FC")
    _bold_to_fc_fcd(
        parc_bold_clean,
        sub_fc_dir,
        dataset,
        parcellation_name,
        config_str,
        save_dynamicFC,
    )


def _bold_to_fc_fcd(
    parc_bold_clean,
    sub_fc_dir,
    dataset,
    parcellation_name,
    config_str="",
    save_dynamicFC=False,
    exc_subcortex=False,
):
    print("calculating FC and FCD")
    os.makedirs(sub_fc_dir, exist_ok=True)
    fc_paths = []
    if exc_subcortex:
        structures_lists = [["ctx"]]
    else:
        structures_lists = [["ctx"], ["sctx", "ctx"]]
    # FC and FCD calculation
    for structures in structures_lists:
        for hemis in [["L"], ["R"], ["L", "R"]]:
            for exc_interhemispheric in [True, False]:
                if exc_interhemispheric & len(hemis) == 1:
                    # drop interhemispheric only makes sense in L-R case
                    continue
                # concatenate parcels (L sctx -> L ctx -> R sctx -> R ctx)
                in_parc_bolds = []
                hemi_parcs = {"L": [], "R": []}  # for exc_interhemispheric
                for hemi in hemis:
                    for structure in structures:
                        hemi_parcs[hemi] += parc_bold_clean[structure][
                            hemi
                        ].index.tolist()
                        in_parc_bolds.append(parc_bold_clean[structure][hemi])
                in_parc_bold = pd.concat(in_parc_bolds, axis=0)
                # FC
                fc = pd.DataFrame(
                    np.corrcoef(in_parc_bold),
                    index=in_parc_bold.index,
                    columns=in_parc_bold.index,
                )
                # excluded interhemispheric connections if indicated
                if exc_interhemispheric:
                    fc.loc[hemi_parcs["L"], hemi_parcs["R"]] = np.NaN
                    fc.loc[hemi_parcs["R"], hemi_parcs["L"]] = np.NaN
                out_path = os.path.join(
                    sub_fc_dir,
                    f'{"_".join(structures)}_parc-{parcellation_name}_hemi-{"".join(hemis)}'
                    + config_str
                    + ("_exc-inter" if exc_interhemispheric else ""),
                )
                fc.to_csv(out_path + "_desc-FC.csv.gz", compression="gzip")
                fc_paths.append(out_path + "_desc-FC.csv.gz")
                fig, ax = plt.subplots(1, figsize=(4, 4))
                sns.heatmap(
                    fc.values, cmap="RdBu_r", vmin=-1, vmax=1, cbar=False, ax=ax
                )
                ax.axis("off")
                fig.tight_layout()
                fig.savefig(out_path + "_desc-FC.png", dpi=80)
                plt.close(fig)
                fc_tril = fc.values[np.tril_indices_from(fc, -1)]
                if exc_interhemispheric:
                    # drop interhemispheric connections (NaNs)
                    fc_tril = fc_tril[~np.isnan(fc_tril)]
                np.savetxt(out_path + "_desc-FCtril.txt", fc_tril, fmt="%.8f")
                # FCD
                fcd_matrix, window_fcs = calculate_fcd(
                    in_parc_bold,
                    window_size=WINDOW_SIZE[dataset],
                    step=STEP[dataset],
                    exc_interhemispheric=exc_interhemispheric,
                    hemi_parcs=hemi_parcs,
                )
                np.savetxt(out_path + "_desc-FCD.txt", fcd_matrix)
                fig, ax = plt.subplots(1, figsize=(4, 4))
                sns.heatmap(fcd_matrix, cbar=False, ax=ax)
                ax.axis("off")
                fig.tight_layout()
                fig.savefig(out_path + "_desc-FCD.png", dpi=80)
                plt.close(fig)
                fcd_tril = fcd_matrix[np.tril_indices_from(fcd_matrix, -1)]
                np.savetxt(out_path + "_desc-FCDtril.txt", fcd_tril, fmt="%.8f")
                if save_dynamicFC:
                    np.savez_compressed(
                        out_path + "_desc-dynamicFC.npz", window_fcs=window_fcs
                    )
    return fc_paths


if __name__ == "__main__":
    dataset = sys.argv[1]
    participant_label = sys.argv[2]
    try:
        session = sys.argv[3]
    except:
        session = None
    sessions = [session]
    parcellations = ["schaefer-100", "schaefer-200"]
    for parc in parcellations:
        for ses in sessions:
            post_fmriprep(
                dataset,
                participant_label,
                ses,
                parcellation_name=parc,
                min_good_scan_duration=MIN_GOOD_SCAN_DURATION[dataset],
            )
