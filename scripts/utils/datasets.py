import os
import numpy as np
import pandas as pd
import scipy.stats
import nibabel
import nilearn.image
import nilearn.surface
import nilearn.signal
from nilearn.input_data import NiftiLabelsMasker
import neuromaps

CODE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
SRC_DIR = os.path.join(CODE_DIR, "src")
from . import transform


def load_pet_map(receptor, parc, zscore=True):
    """
    Preprocesses PET maps by Z-scoring them and taking a
    weighted average in case multiple maps exist for a given
    receptor x tracer combination

    Parameters
    ----------
    receptor: {"nmda", "gabaa"}
    parc: (str)
    zscore: (bool)

    Returns
    -------
    pet_parcellated: (pd.Series)
    """
    pet_filenames = {
        "nmda": "NMDA_ge179_hc29_galovic.nii.gz",
        "gabaa": "GABAa-bz_flumazenil_hc16_norgaard.nii.gz",
    }
    pet_img = os.path.join(SRC_DIR, "PET_nifti_images", pet_filenames[receptor])
    masker = NiftiLabelsMasker(
        os.path.join(
            SRC_DIR, "parcellations", f"tpl-MNI152_desc-{parc}_parcellation_1mm.nii.gz"
        ),
        strategy="sum",
        resampling_target="data",
        background_label=0,
    )
    # count the number of non-zero voxels per parcel so the average
    # is calculated only among non-zero voxels (visualizing the PET map
    # on volumetric parcellations, the parcels are usually much thicker
    # than the PET map on the cortex, and there are a large number of
    # zero PET values in each parcel which can bias the parcelled values)
    nonzero_mask = nilearn.image.math_img("pet_img != 0", pet_img=pet_img)
    nonzero_voxels_count_per_parcel = masker.fit_transform(nonzero_mask).flatten()
    # take the average of PET values across non-zero voxels
    pet_value_sum_per_parcel = masker.fit_transform(pet_img).flatten()
    pet_parcellated = pet_value_sum_per_parcel / nonzero_voxels_count_per_parcel
    pet_parcellated = pd.Series(
        pet_parcellated, index=transform.load_ordered_parcel_labels(parc)
    )
    if zscore:
        pet_parcellated = scipy.stats.zscore(pet_parcellated)
    return pet_parcellated


def load_maps(parcellation_name, maps="6maps", zscore=True):
    """
    Loads the heterogeneity maps

    Parameters
    ---------
    parcellation_name: (str)
    maps: (str or list of str)
        - '6maps': the six maps used in the models
        - 'sydnor2021': maps shown in Fig2 of Sydnor 2021 (https://doi.org/10.1016/j.neuron.2021.06.016)
         except externopyramidization is replaced by LTC G1 and for evoexp xu2020 is used
        - '6maps_sydnor2021': combination of the above
        - list: a list including any valid neuromaps / receptor / other maps
    zscore: (bool)
    """
    if maps == "6maps":
        maps = ["myelinmap", "thickness", "fcgradient01", "genepc1", "nmda", "gabaa"]
    elif maps == "sydnor2021":
        maps = [
            "SAaxis",
            "myelinmap",
            "fcgradient01",
            "evoexp",
            "scalingpnc",
            "cmr02",
            "cmruglu",
            "meancbf",
            "genepc1",
            "cogpc1",
            "ltcg1",
            "thickness",
        ]
    elif maps == "6maps_sydnor2021":
        maps = [
            "myelinmap",
            "thickness",
            "fcgradient01",
            "genepc1",
            "nmda",
            "gabaa",
            "SAaxis",
            "evoexp",
            "scalingpnc",
            "cmr02",
            "cmruglu",
            "meancbf",
            "cogpc1",
            "ltcg1",
        ]
    maps_parcellated = []
    for curr_map in maps:
        if curr_map in ["nmda", "gabaa"]:
            maps_parcellated.append(
                load_pet_map(curr_map, parcellation_name, zscore=zscore).rename(
                    curr_map
                )
            )
        else:
            if curr_map == "ltcg1":
                map_paths = [
                    os.path.join(
                        SRC_DIR, "tpl-fs_LR_hemi-L_den-32k_desc-LTC_G1.shape.gii"
                    ),
                    os.path.join(
                        SRC_DIR, "tpl-fs_LR_hemi-R_den-32k_desc-LTC_G1.shape.gii"
                    ),
                ]
                space = "fsLR"
                den = "32k"
            else:
                source = None
                if curr_map == "evoexp":
                    source = "xu2020"
                # get space
                annotations = neuromaps.datasets.available_annotations(
                    desc=curr_map, source=source
                )
                if len(annotations) > 1:
                    raise NotImplementedError
                space = annotations[0][2]
                den = annotations[0][3]
                # get data
                map_paths = neuromaps.datasets.fetch_annotation(
                    desc=curr_map, source=source
                )
                # make sure L is before R
                if isinstance(map_paths, list) and (len(map_paths) > 1):
                    map_paths = sorted(map_paths)
                    assert "hemi-L" in map_paths[0]
            if space == "MNI152":
                parcellated = transform.parcellate_vol(
                    map_paths, parcellation_name, concat=True
                ).iloc[:, 0]
            else:
                # for fsLR if it is 164k downsample it as the parcellation maps are
                # available in 32k density
                if (space == "fsLR") & (den != "32k"):
                    maps_32k = neuromaps.transforms.fslr_to_fslr(map_paths, "32k")
                    map_data = {
                        "L": maps_32k[0].agg_data(),
                        "R": maps_32k[1].agg_data(),
                    }
                    den = "32k"
                # transform civet to fsLR
                elif space == "civet":
                    maps_32k = neuromaps.transforms.civet_to_fslr(map_paths, "32k")
                    map_data = {
                        "L": maps_32k[0].agg_data(),
                        "R": maps_32k[1].agg_data(),
                    }
                    space = "fsLR"
                    den = "32k"
                else:
                    map_data = {
                        "L": nibabel.load(map_paths[0]).agg_data(),
                        "R": nibabel.load(map_paths[1]).agg_data(),
                    }
                # determine downsampling (fsaverage 10k or fsLR 32k)
                if space == "fsaverage":
                    downsampled = den == "10k"
                elif space == "fsLR":
                    downsampled = den == "32k"
                # parcellate
                parcellated = transform.parcellate_surf(
                    map_data,
                    parcellation_name,
                    space=space,
                    downsampled=downsampled,
                    concat=True,
                ).iloc[:, 0]
            if zscore:
                parcellated = scipy.stats.zscore(parcellated)
            annotations = neuromaps.datasets.available_annotations(
                desc=curr_map, source=source
            )
            if len(annotations) > 1:
                raise NotImplementedError
            maps_parcellated.append(parcellated.rename(curr_map))
    maps_parcellated = pd.concat(maps_parcellated, axis=1)
    return maps_parcellated


def load_mesh_paths(kind="orig", space="fsaverage", downsampled=True):
    """
    Loads or creates surfaces of left and right hemispheres

    Parameters
    ----------
    kind: (str)
        - 'orig'
        - 'inflated'
        - 'sphere'
    space: (str)
        - 'fsaverage'
        - 'fs_LR'
    downsampled: (bool)

    Returns
    ------
    paths: (dict of str) path to downsampled surfaces of L and R
    """
    paths = {}
    for hem in ["L", "R"]:
        if space == "fsaverage":
            hem_fullname = {"L": "left", "R": "right"}
            if downsampled:
                fsa_version = "fsaverage5"
            else:
                fsa_version = "fsaverage"
            if kind == "orig":
                paths[hem] = nilearn.datasets.fetch_surf_fsaverage(fsa_version)[
                    f"pial_{hem_fullname[hem]}"
                ]
            elif kind == "inflated":
                paths[hem] = nilearn.datasets.fetch_surf_fsaverage(fsa_version)[
                    f"infl_{hem_fullname[hem]}"
                ]
            elif kind == "sphere":
                paths[hem] = nilearn.datasets.fetch_surf_fsaverage(fsa_version)[
                    f"sphere_{hem_fullname[hem]}"
                ]
            elif kind == "semi-inflated":
                if fsa_version == "fsaverage5":
                    raise ValueError(
                        "Semi-inflated surfaces are not available for fsaverage5"
                    )
                else:
                    paths[hem] = os.path.join(
                        SRC_DIR, f"fsaverage_{hem.lower()}h.pial_semi_inflated"
                    )
        elif space == "fs_LR":
            if kind == "orig":
                paths[hem] = os.path.join(
                    SRC_DIR, f"S1200.{hem}.pial_MSMAll.32k_fs_LR.surf.gii"
                )
            else:
                raise NotImplementedError
    return paths


def load_ahba_data(parcellation_name):
    """
    Loads parcellated AHBA data previously processed
    using `abagen` and https://github.com/amnsbr/laminar_organization
    scripts
    """
    assert parcellation_name == "schaefer-100", "Only schaefer-100 is supported"
    file_path = os.path.join(
        SRC_DIR, "ahba_parc-schaefer100_hemi-L_ibf_threshold-0.5_missing-centroids.npz"
    )
    return np.load(file_path, allow_pickle=True)["data"].tolist()["all"]


def load_aggregate_gene_expression(gene_list, parcellation_name, avg_method="mean"):
    """
    Gets the aggregate expression of genes in `gene_list`

    Parameters
    ---------
    gene_list: (list | pd.Series)
        list of gene names or a series with names in the index
        and weights in the values
        Note: This will ignore genes that do not exist in the current
        version of gene expression data
    parcellation_name: (str)
    avg_method: (str)
        - mean
        - median (ignores the weights in gene list)
    """
    # get the ahba data
    ahba_data = load_ahba_data(parcellation_name)
    # get the gene list that exist in ahba data
    ahba_genes = list(ahba_data.columns)
    if isinstance(gene_list, list):
        # if no weights are given set the weight
        # of all genes to 1
        gene_list = pd.Series(1, index=gene_list)
    exist_gene_list = set(gene_list.index) & set(ahba_genes)
    print(
        f"{gene_list.shape[0] - len(exist_gene_list)} of {gene_list.shape[0]} genes do not exist"
    )
    gene_list = gene_list.loc[exist_gene_list]
    # get the aggregate expression
    if avg_method == "mean":
        aggregate_expression = (
            ahba_data.loc[:, gene_list.index]  # gene expression
            @ gene_list.values  # weights
        ) / gene_list.sum()
    elif avg_method == "median":
        aggregate_expression = ahba_data.loc[:, gene_list.index].median(axis=1)
    return aggregate_expression
