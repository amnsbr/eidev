import os
import numpy as np
import pandas as pd
import nibabel
from nilearn.input_data import NiftiLabelsMasker

CODE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
SRC_DIR = os.path.join(CODE_DIR, "src")

N_VERTICES_HEM_FSLR = 32492
N_VERTICES_HEM_ICO5 = 10242
N_VERTICES_HEM_ICO7 = 163842

MIDLINE_PARCELS = {
    "schaefer-100": ["Background+FreeSurfer_Defined_Medial_Wall"],
    "schaefer-200": ["Background+FreeSurfer_Defined_Medial_Wall"],
}

def load_parcellation_map(
    parcellation_name,
    concatenate=False,
    load_indices=False,
    space="fsaverage5",
    downsampled=False,
):
    """
    Loads parcellation maps of L and R hemispheres in fsa5 space, correctly relabels them
    and concatenates them

    Parameters
    ----------
    parcellation_name: (str)
    concatenate: (bool)
        cocnatenate the hemispheres
    load_indices: (bool)
        return the parcel indices instead of their names
    space: {fsaverage, fsaverage5, fsLR}
    downsampled: (bool)
        with 'fsaverage', if true, will change the space to 'fsaverage5'
        this is just for compatiability with stats.get_rotated_parcels

    Returns
    -------
    parcellation_map: (np.ndarray or dict of np.ndarray)
    """
    # fsLR
    if space == "fsLR":
        n_parcels = int(parcellation_name.split("-")[1])
        parcellation_path = os.path.join(
            SRC_DIR,
            "parcellations",
            f"Schaefer2018_{n_parcels}Parcels_7Networks_order.dlabel.nii",
        )
        parcellation_map = nibabel.load(parcellation_path).get_fdata().squeeze()
        if not load_indices:
            labels = load_ordered_parcel_labels(parcellation_name).tolist()
            labels = ["Background+FreeSurfer_Defined_Medial_Wall"] + labels
            transdict = dict(enumerate(labels))
            parcellation_map = np.vectorize(transdict.get)(parcellation_map)
        if concatenate:
            return parcellation_map
        else:
            return {
                "L": parcellation_map[:N_VERTICES_HEM_FSLR],
                "R": parcellation_map[N_VERTICES_HEM_FSLR:],
            }
    # fsaverage / fsaverage5
    if (space == "fsaverage") & (downsampled):
        space = "fsaverage5"
    parcellation_map = {}
    for hem in ["L", "R"]:
        if space == "fsaverage5":
            parcellation_path = os.path.join(
                SRC_DIR, "parcellations", f"{hem.lower()}h.{parcellation_name}.annot"
            )
        elif space == "fsaverage":
            parcellation_path = os.path.join(
                SRC_DIR,
                "parcellations",
                f"{hem.lower()}h.{parcellation_name}.fsa.annot",
            )
        parcellation_map[hem], _, labels = nibabel.freesurfer.io.read_annot(
            parcellation_path
        )
        # label parcellation map if indicated
        if not load_indices:
            # labels post-processing for each specific parcellation
            labels = list(map(lambda l: l.decode(), labels))  # b'name' => 'name'
            transdict = dict(enumerate(labels))
            parcellation_map[hem] = np.vectorize(transdict.get)(parcellation_map[hem])
    if concatenate:
        return np.concatenate([parcellation_map["L"], parcellation_map["R"]])
    else:
        return parcellation_map


def load_ordered_parcel_labels(parcellation_name):
    """
    Loads the lables corresponding to the parcels in volumetric
    space, except for the parcel ID 0, which corresponds to
    background / midline

    Parameter
    --------
    parcellation_name: (str)

    Returns
    -------
    labels: (np.array)
    """
    n_parcels = int(parcellation_name.replace("schaefer-", ""))
    lut_path = os.path.join(
        SRC_DIR, "parcellations", f"Schaefer2018_{n_parcels}Parcels_7Networks_order.txt"
    )
    labels = pd.read_csv(lut_path, sep="\t", header=None)[1].values
    return labels


def parcellate_surf(
    surface_data,
    parcellation_name,
    method="mean",
    midline="drop",
    space="fsaverage5",
    align_order=True,
    concat=False,
    downsampled=False,
):
    """
    Parcellates `surface data` using `parcellation` and by taking the
    median or mean (specified via `averaging_method`) of the vertices within each parcel.

    Parameters
    ----------
    surface_data: (dict of np.ndarray)
        n_vertices x n_features surface data of L and R hemispheres
    parcellation_name: (str)
    method: (str)
        Method of aggregating over vertices within a parcel.
            - 'mean' [default]
            - 'median'
            - 'sum'
    midline: (str | None)
        - None: keeps midline as is
        - drop: drops midline [default]
        - nan: sets midline to nan
    space: (str)
        used if dataset & participant_label are None
        - 'fsaverage5' [default]
        - 'fsaverage'
        - 'fsLR'
    align_order: (bool)
        align parcel order to other parcellated data (based on micapipe lut files)
    concat: (bool)
    downsampled: (bool)
        corresponds to 10k in fsaverage (fsaverage5) and 32k in fsLR

    Returns
    ---------
    parcellated_data: (pd.DataFrame or dict of pd.DataFrame) n_parcels x n_features for data of L and R hemispheres or both hemispheres
    """
    # load parcellation map
    labeled_parcellation_maps = load_parcellation_map(
        parcellation_name, concatenate=False, space=space, 
        downsampled=downsampled
    )
    # load micapipe lut for correct ordering of parcels
    lut = pd.read_csv(
        os.path.join(SRC_DIR, "parcellations", f"lut_{parcellation_name}_mics.csv")
    ).set_index("label")
    lut_cortical = lut.loc[(lut["mics"] >= 1000)]
    parcellated_data = {}
    for hem in ["L", "R"]:
        # parcellate
        parcellated_vertices = pd.DataFrame(
            surface_data[hem], index=labeled_parcellation_maps[hem]
        )
        parcellated_vertices = parcellated_vertices.reset_index(drop=False).groupby(
            "index"
        )
        # operate on groupby object if needed
        if method == "median":
            parcellated_data[hem] = parcellated_vertices.median()
        elif method == "mean":
            parcellated_data[hem] = parcellated_vertices.mean()
        elif method == "sum":
            parcellated_data[hem] = parcellated_vertices.sum()
        # remove midline data
        if midline == "drop":
            parcellated_data[hem] = parcellated_data[hem].drop(
                index=parcellated_data[hem].index.intersection(
                    MIDLINE_PARCELS[parcellation_name]
                )
            )
        elif midline == "nan":
            parcellated_data[hem].loc[
                parcellated_data[hem].index.intersection(
                    MIDLINE_PARCELS[parcellation_name]
                )
            ] = np.NaN
        if align_order:
            # correctly order parcels
            ordered_parcels = lut_cortical.index.intersection(
                parcellated_data[hem].index
            )
        else:
            ordered_parcels = parcellated_data[hem].index
        parcellated_data[hem] = parcellated_data[hem].loc[ordered_parcels]
    if concat:
        return pd.concat([parcellated_data["L"], parcellated_data["R"]])
    else:
        return parcellated_data


def deparcellate_surf(
    parcellated_data, parcellation_name, space="fsaverage5", concat=False
):
    """
    Project the parcellated data to surface vertices while handling empty parcels
    (parcels that are not in the parcellated data but are in the parcellation map)

    Parameters
    ----------
    parcellated_data: (pd.DataFrame | pd.Series) n_parcels x n_features
    parcellation_name: (str)
    dataset: (str)
    participant_label: (str)
    ses: (str)
    space: (str)
        used when dataset & participant_label are None
        - 'fsaverage5' [default]
        - 'fsaverage'
        - 'fsLR'

    Returns
    -------
    surface_map: (np.ndarray) n_vertices [both hemispheres] x n_features
    """
    # load concatenated parcellation map
    parcellation_map = load_parcellation_map(
        parcellation_name, concatenate=False, space=space
    )
    surface_data = {}
    # load dummy parcellated data covering the whole brain
    dummy_surf_data = {
        "L": np.zeros_like(parcellation_map["L"]),
        "R": np.zeros_like(parcellation_map["R"]),
    }
    parcellated_dummy = parcellate_surf(
        dummy_surf_data, parcellation_name, midline=None, align_order=False, space=space
    )
    for hemi in ["L", "R"]:
        all_parcels = parcellated_dummy[hemi].index.to_series().rename("parcel")
        # create a dataframe including all parcels, where invalid parcels are NaN
        #   (this is necessary to be able to project it to the parcellation)
        labeled_parcellated_data = pd.concat(
            [parcellated_data, all_parcels], axis=1
        ).set_index("parcel")
        # get the surface map by indexing the parcellated map at parcellation labels
        surface_data[hemi] = labeled_parcellated_data.loc[
            parcellation_map[hemi]
        ].values.squeeze()
    if concat:
        surface_data = np.concatenate([surface_data["L"], surface_data["R"]])
    return surface_data


def parcellate_vol(
    img,
    parcellation_name,
    method="mean",
    midline="drop",
    align_order=True,
    concat=False,
):
    """
    Parcellated volumetric image. Note that volumetric image must be in MNI152 space.
    """
    masker = NiftiLabelsMasker(
        os.path.join(
            SRC_DIR,
            "parcellations",
            f"tpl-MNI152_desc-{parcellation_name}_parcellation_1mm.nii.gz",
        ),
        strategy=method,
        resampling_target="data",
        background_label=0,
    )
    parcellated_data = masker.fit_transform(img)
    parcellated_data = pd.DataFrame(
        parcellated_data.T, index=load_ordered_parcel_labels(parcellation_name)
    )
    # drop or NaN midline
    if midline == "drop":
        parcellated_data = parcellated_data.drop(
            index=parcellated_data.index.intersection(
                MIDLINE_PARCELS[parcellation_name]
            )
        )
    elif midline == "nan":
        parcellated_data.loc[
            parcellated_data.index.intersection(MIDLINE_PARCELS[parcellation_name])
        ] = np.NaN
    # align order with mica lut
    if align_order:
        lut = pd.read_csv(
            os.path.join(
                CODE_DIR, "src", "parcellations", f"lut_{parcellation_name}_mics.csv"
            )
        ).set_index("label")
        lut_cortical = lut.loc[(lut["mics"] >= 1000)]
        ordered_parcels = lut_cortical.index.intersection(parcellated_data.index)
    else:
        ordered_parcels = parcellated_data.index
    parcellated_data = parcellated_data.loc[ordered_parcels]
    # break into hemis if needed
    if concat:
        return parcellated_data
    else:
        return {
            "L": parcellated_data.iloc[: parcellated_data.shape[0] // 2],
            "R": parcellated_data.iloc[parcellated_data.shape[0] // 2 :],
        }
