import os
import time
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats
import requests
import nilearn.surface
import nilearn.plotting
import enigmatoolbox.permutation_testing
from sklearn.cross_decomposition import PLSRegression
import statsmodels
from tqdm import tqdm

from . import datasets, transform

DSEA_STAGES = [
    "Early Fetal",
    "Early Mid Fetal",
    "Late Mid Fetal",
    "Late Fetal",
    "Neotal Early Infancy",
    "Late Infancy",
    "Early Childhood",
    "Middle Late Childhood",
    "Adolescence",
    "Young Adulthood",
]


def get_rotated_parcels(
    parcellation_name,
    n_perm,
    excluded_parcels=[],
    return_indices=True,
    space="fsaverage",
    downsampled=True,
    seed=0,
):
    """
    Uses ENIGMA Toolbox approach to spin parcel centroids on the cortical sphere instead
    of spinning all the vertices.

    Parameters
    ----------
    parcellation_name: (str)
    n_perm: (int)
    excluded_parcels: (list)
    return_indices: (bool)
        returns rotated indices instead of labels
    space: {'fsaverage'}
    downsampled: (bool)
    seed: (int)
    """
    # get the coordinates for the centroids of each parcel
    # on cortical sphere of the space
    centroids = {}
    for hem in ["L", "R"]:
        surf_path = datasets.load_mesh_paths(
            "sphere", space=space, downsampled=downsampled
        )[hem]
        vertices = nilearn.surface.load_surf_mesh(surf_path).coordinates
        parc = transform.load_parcellation_map(
            parcellation_name, False, downsampled=downsampled, space=space
        )[hem]
        centroids[hem] = pd.DataFrame(columns=["x", "y", "z"])
        for parcel_name in np.unique(parc):
            if parcel_name not in (
                excluded_parcels + transform.MIDLINE_PARCELS.get(parcellation_name, [])
            ):
                this_parc = np.where(parc == parcel_name)[0]
                centroids[hem].loc[parcel_name] = np.mean(
                    vertices[this_parc, :], axis=0
                )
    # rotate the parcels
    np.random.seed(seed)
    rotated_indices = enigmatoolbox.permutation_testing.rotate_parcellation(
        centroids["L"].values.astype("float"),
        centroids["R"].values.astype("float"),
        n_perm,
    ).astype("int")
    # get the rotated parcel names
    rotated_parcels = pd.concat(
        [centroids["L"], centroids["R"]], axis=0
    ).index.to_numpy()[rotated_indices]
    if return_indices:
        # current indices are based on parcels that are L->R
        # and then alphabetically ordered, which may be different
        # from the original order (e.g. in schaefer). Map
        # the parcel names to their indices in correctly-ordered
        # parcellation map
        ordered_parcels = enumerate(
            transform.load_ordered_parcel_labels(parcellation_name)
        )
        ordered_parcels = {parcel: i for i, parcel in ordered_parcels}
        rotated_indices = np.vectorize(ordered_parcels.get)(rotated_parcels)
        return rotated_indices
    else:
        return rotated_parcels


def spin_test_parcellated(
    X, Y, parcellation_name, n_perm=1000, space="fsaverage", seed=0, method="pearson",
):
    """
    Uses ENIGMA Toolbox approach to spin parcel centroids on the cortical sphere instead
    of spinning all the vertices. Ideally spin the parcellated data with less NaN values
    by assigining it to X.

    Parameters
    ---------
    X, Y: (pd.DataFrame) n_parcel x n_features
    parcellation_name: (str)
    n_perm: (int)
    space: {'fsaverage'}
    seed: (int)
    method: {'pearson', 'kendall', 'spearman'} or callable
        passed on to `pandas.DataFrame.corr`
    """
    # calculate test correlation coefficient between all pairs of columns between surface_data_to_spin and surface_data_target
    coefs = (
        pd.concat([X, Y], axis=1)
        # calculate the correlation coefficient between all pairs of columns within and between X and Y
        .corr(method=method)
        # select only the correlations we are interested in
        .iloc[: X.shape[1], -Y.shape[1] :]
        # convert it to shape (1, n_features_Y, n_features_surface_X)
        .T.values[np.newaxis, :]
    )
    # get the spin rotated parcel indices
    rotated_indices = get_rotated_parcels(
        parcellation_name, n_perm, return_indices=True, space=space, seed=seed
    )
    # add NaN parcels back to X and Y so that the number of parcels
    # in them and rotated parcels is the same (to make sure the unlabeled
    # numpy arrays are aligned)
    X_all_parcels = transform.parcellate_surf(
        transform.deparcellate_surf(X, parcellation_name, space=space, concat=False),
        parcellation_name,
        space=space,
        concat=True,
    )
    X_all_parcels = X_all_parcels.loc[
        ~X_all_parcels.index.isin(transform.MIDLINE_PARCELS.get(parcellation_name, []))
    ]
    X_all_parcels.columns = X.columns
    Y_all_parcels = transform.parcellate_surf(
        transform.deparcellate_surf(Y, parcellation_name, space=space, concat=False),
        parcellation_name,
        space=space,
        concat=True,
    )
    Y_all_parcels = Y_all_parcels.loc[
        ~Y_all_parcels.index.isin(transform.MIDLINE_PARCELS.get(parcellation_name, []))
    ]
    Y_all_parcels.columns = Y.columns
    assert X_all_parcels.shape[0] == Y_all_parcels.shape[0] == rotated_indices.shape[0]
    null_distribution = np.zeros((n_perm, Y.shape[1], X.shape[1]))
    for x_col in range(X_all_parcels.shape[1]):
        # get all surrogate parcellated maps at once.
        # this involves complicated indexing but the is best way
        # to achieve this extremely more efficiently than loops.
        surrogates = np.zeros((X_all_parcels.shape[0], n_perm))
        surrogates[
            rotated_indices.T.flatten(),
            np.arange(0, n_perm)
            .reshape(n_perm, 1)
            .repeat(repeats=X_all_parcels.shape[0], axis=1)
            .flatten(),
        ] = np.tile(X_all_parcels.values[:, x_col], n_perm)
        surrogates = pd.DataFrame(surrogates)
        surrogates.columns = [f"surrogate_{i}" for i in range(n_perm)]
        null_distribution[:, :, x_col] = (
            pd.concat([surrogates, Y_all_parcels.reset_index(drop=True)], axis=1)
            .corr(method=method)
            .iloc[: surrogates.shape[1], -Y.shape[1] :]
            .values
        )
    # calculate p value
    pvals = (np.abs(null_distribution) > np.abs(coefs)).mean(axis=0)
    # remove unnecessary dimension of test_r
    coefs = coefs[0, :, :]
    pvals = pd.DataFrame(pvals, index=Y.columns, columns=X.columns)
    coefs = pd.DataFrame(coefs, index=Y.columns, columns=X.columns)
    return coefs, pvals, null_distribution

def anova(parc_data, categories, output="text", force_posthocs=False):
    """
    Compares value of parcellated data across `categories`
    using ANOVA and post-hoc t-tests

    Parameters
    ----------
    parc_data: (pd.DataFrame) (n_parcel, n_features)
    categories: (pd.Series or np.ndarray) (n_parcel,)
    output: {'text', 'stats'}
    force_posthocs: (bool)
        forces posthocs to be done regardless of
        ANOVA p value

    Returns
    -------
    anova_res_str: (str)
    OR
    anova_res: (pd.Series)
    """
    F, p_val = scipy.stats.f_oneway(
        *[
            category_data[1].dropna().values
            for category_data in parc_data.groupby(categories)
        ]
    )
    anova_res_str = f"----\nF statistic {F}, pvalue {p_val}\n"
    anova_res = pd.Series({"F": F})
    if (p_val < 0.05) or force_posthocs:
        # set alpha for bonferroni correction across all pairs of categories
        alpha = 0.05 / len(list(itertools.combinations(categories.cat.categories, 2)))
        if force_posthocs:
            anova_res_str += f"\tPost-hoc T-tests:\n(bonferroni alpha: {alpha})\n"
        else:
            anova_res_str += f"\tPost-hoc T-tests passing alpha of {alpha}:\n"
        for cat1, cat2 in itertools.combinations(categories.cat.categories, 2):
            t_statistic, t_p = scipy.stats.ttest_ind(
                parc_data.loc[categories == cat1].dropna(),
                parc_data.loc[categories == cat2].dropna(),
            )
            anova_res.loc[f"{cat1}-{cat2}"] = t_statistic
            if (t_p < alpha) or (force_posthocs):
                anova_res_str += f"\t\t{cat1} - {cat2}: T {t_statistic}, p {t_p}\n"
    if output == "text":
        return anova_res_str
    else:
        return anova_res

def anova_spin(parc_data, categories, parcellation_name, n_perm=1000, seed=921):
    """
    Performs anova and posthoc t-tests using spin permutation

    Parameters
    ----------
    parc_data: (pd.DataFrame) (n_parcel, n_features)
    categories: (pd.Series or np.ndarray) (n_parcel,)
    parcellation_name: (str)
    n_perm: (int)
    seed: (int)
    """
    # get the test statistics
    test_stats = anova(parc_data, categories, output="stats", force_posthocs=True)
    # create spin surrogates
    rotated_parcels = get_rotated_parcels(
        parcellation_name, n_perm, return_indices=False, seed=seed,
    )
    # create null distribution
    null_dist = np.zeros((test_stats.shape[0], n_perm))
    for i in tqdm(range(n_perm)):
        surrogate_data = parc_data.loc[rotated_parcels[:, i]]
        surrogate_data.index = parc_data.index
        null_dist[:, i] = anova(
            surrogate_data, categories, output="stats", force_posthocs=True
        ).values
    # p-value
    p_vals = (np.abs(null_dist) > np.abs(test_stats[:, np.newaxis])).mean(axis=1)
    p_vals = pd.Series(p_vals, index=test_stats.index)
    return test_stats, p_vals

def madicc(x, y):
    """
    Median Absolute Deviation Intraclass Correlation Coefficient

    This function implements the intraclass version of the Median Absolute
    Deviation Correlation Coefficient as described in Shevlyakov & Smirnov (2011).

    Translated from the original Matlab code in the following link:
    https://warwick.ac.uk/fac/sci/statistics/staff/academic-research/nichols/scripts/matlab/madicc.m

    Parameters
    ----------
    x, y: (np.ndarray) of shape (n,)

    Returns
    -------
    rmad: float
    """
    I = ~np.isnan(x) & ~np.isnan(y)
    if not I.any():
        rmad = np.nan
    else:
        mx = np.median(x[I])
        my = np.median(y[I])
        Sp = (x[I] - mx) + (y[I] - my)
        Sm = (x[I] - mx) - (y[I] - my)
        madSp = np.median(np.abs(Sp - np.median(Sp)))
        madSm = np.median(np.abs(Sm - np.median(Sm)))
        if madSp == 0 and madSm == 0:
            rmad = np.nan
        else:
            rmad = (madSp**2 - madSm**2) / (madSp**2 + madSm**2)
    return rmad

def ahba_pls(parcellated_data, parcellation_name='schaefer-100', n_genes=500, n_components=1):
    """
    Performs PLS between columns of `parc_data` and AHBA gene expression data

    Parameters
    ---------
    parcellated_data: (pd.DataFrame)
    parcellation_name: (str)
    n_genes: (int)
        number of genes with the highest absolute weight to select
    n_components: (int)
        number of PLS components

    Returns
    -------
    top_genes: (list of dict)
        each element corresponds to a PLS component and includes
        list of top genes with 'pos' and 'neg' weights
    pls: (sklearn.cross_decomposition.PLSRegression)
        fitted pls model
    """
    ahba_df = datasets.load_ahba_data(parcellation_name)
    shared_parcels = ahba_df.index.intersection(parcellated_data.index)
    x = ahba_df.loc[shared_parcels]
    y = parcellated_data.loc[shared_parcels]
    pls = PLSRegression(n_components=n_components)
    pls.fit(x, y)
    weights=pd.DataFrame(pls.x_weights_, index=ahba_df.columns)
    top_genes = []
    for i in range(n_components):
        curr_top_genes = weights.iloc[:, i].abs().sort_values(ascending=False).iloc[:n_genes].index.tolist()
        weights_top_genes = weights.iloc[:, i].loc[curr_top_genes]
        top_pos = weights_top_genes.loc[weights_top_genes >= 0].index.tolist()
        top_neg = weights_top_genes.loc[weights_top_genes < 0].index.tolist()
        top_genes.append({'pos': top_pos, 'neg': top_neg})
    return top_genes, pls

def fetch_dsea_results(gene_list, fdr=True, mirror='new'):
    """
    Runs CSEA tool developmental enrichment on the gene list 
    using http://doughertylab.wustl.edu/csea-tool-2/, parses
    and cleans the results table and returns it

    Parameters
    ---------
    gene_list: (list of str)
    fdr: (bool)
        fetch fdr-corrected p-values
    mirror: (str)
        online app url

    Returns
    -------
    csea_res: (pd.DataFrame)
        CSEA developmental enrichment output table
    """
    if mirror == 'new':
        url = 'http://doughertylab.wustl.edu/jdlab-fisher/cgi-bin/seaBrainRegion.cgi'
    elif mirror == 'old':
        url = 'http://doughertytools.wustl.edu/cgi-bin/seaBrainRegion.cgi'
    else:
        url = mirror
    res = requests.post(url, 
                {
                    'symbols': ' '.join(gene_list), 
                    'Candidate gene list from': 'Human'
                })
    if res.status_code != 200:
        print("Request failed. Try another mirror")
        return
    csea_res = pd.read_html(res.text)[0]
    csea_res.columns = csea_res.iloc[0, :]
    csea_res = csea_res.drop(0).set_index('Brain Regions and Development and P-Values')
    if fdr:
        csea_res = csea_res.applymap(lambda c: (c.split('(')[-1][:-1])).astype('float')
    else:
        csea_res = csea_res.applymap(lambda c: (c.split('(')[0])).astype('float')
    return csea_res

def dsea(
    X, parcellation_name, seed=0, n_genes=500,
    mirror='new', pSI='0.05', fdr=True
):
    """
    Performs developmentla specific expression analysis (dSEA) on the
    provided parcellated map

    Parameters
    ---------
    X: (pd.DataFrame)
    parcellation_name: (str)
    seed: (int)
        Seed used for PLS
    n_genes: (int)
        Number of top genes to select
    mirror: {'new', 'old', str}
        CSEA tool url
    pSI: specificity index threshold
    fdr: (bool)
        apply FDR correction on p-values

    Returns
    -------
    nlog_dsea_res: (pd.Series)
        negative log10 of p-values
    raw_dsea_res: (dict)
        raw results
    top_genes: (list of dict)
        top genes from PLS with 'pos' and 'neg' weights
    pls: (sklearn.cross_decomposition.PLSRegression)
        fitted PLS model
    """
    # run PLS and get the top associated genes
    np.random.seed(seed)
    top_genes, pls = ahba_pls(X, parcellation_name=parcellation_name, n_genes=n_genes, n_components=1)
    # make sure 'pos' genes correspond to genes with *higher* expression
    # towards regions with *higher* values in X
    # therefore when y weight is negative switch positive and negative
    if pls.y_weights_[0] < 0:
        top_genes[0] = {
            "pos": top_genes[0]["neg"],
            "neg": top_genes[0]["pos"],
        }
    # run online CSEA tool on the positive and negative lists of genes
    pos_res = fetch_dsea_results(top_genes[0]["pos"], fdr=False, mirror=mirror)
    time.sleep(2)  # to avoid back-to-back requests to the online tool
    neg_res = fetch_dsea_results(top_genes[0]["neg"], fdr=False, mirror=mirror)
    raw_dsea_res = {"Positive": pos_res, "Negative": neg_res}
    # filter to cortex, apply FDR on cortex p-values and take negative log10
    nlog_dsea_res = {}
    fdr_sigs = {}
    for k, dsea_res in raw_dsea_res.items():
        dsea_res = dsea_res.copy()
        # select pSI threshold
        dsea_res = dsea_res.loc[:, [pSI]]
        # clean structure and stages names
        dsea_res["Structure"] = dsea_res.index.to_series().apply(
            lambda s: s.split(".")[0].strip()
        )
        dsea_res["Stage"] = dsea_res.index.to_series().apply(
            lambda s: " ".join(s.split(".")[1:])
        )
        dsea_res["Stage"] = dsea_res["Stage"].str.strip()
        dsea_res = (
            dsea_res.reset_index(drop=True)
            .set_index(["Stage", "Structure"])["0.05"]
            .unstack()
        )
        # select cortex and order stages
        dsea_res = dsea_res.loc[DSEA_STAGES, "Cortex"]
        # FDR if indicated
        if fdr:
            fdr_sig, p_fdr = statsmodels.stats.multitest.fdrcorrection(dsea_res.values)
            dsea_res.iloc[:] = p_fdr
            fdr_sigs[k] = fdr_sig
        # take negative log10 of p-values
        dsea_res = -np.log10(dsea_res)
        nlog_dsea_res[k] = dsea_res
    # convert the results to a stacked series
    nlog_dsea_res = pd.DataFrame(nlog_dsea_res).stack()
    return nlog_dsea_res, raw_dsea_res, top_genes, pls

def dsea_plot(
    nlog_dsea_res, sigs={}, ax=None, 
    colors={}, offsets={}, legend_kwargs={}
):
    """
    Plots the dSEA results

    Parameters
    ---------
    nlog_dsea_res: (pd.Series)
    sigs: (dict)
        significant indicators of 'pos' and 'neg' sets
    ax: (matplotlib.axes.Axes)
    colors: (dict)
        colors for 'Positive' and 'Negative' sets
    offsets: (dict)
        offsets for significance indicators
    legend_kwargs: (dict)
        kwargs for legend

    Returns
    -------
    ax: (matplotlib.axes.Axes)
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8.5, 2))
    plot_data = nlog_dsea_res.unstack().unstack().reset_index()
    sns.barplot(
        data=plot_data,
        x="Stage",
        y=0,
        hue="level_0",
        palette=[
            colors.get('Positive', 'red'), 
            colors.get('Negative', 'blue')
        ],
        saturation=100,
        alpha=1,
        ax=ax,
    )
    # set default offsets for significance indicators
    offsets['Positive'] = offsets.get('Positive', -0.2)
    offsets['Negative'] = offsets.get('Negative', +0.2)
    # significance indicators
    for k in ['Positive', 'Negative']:
        curr_nlog_dsea_res = nlog_dsea_res.loc[slice(None), k]
        # get provided significant indicators (e.g.
        # from spin test) and when not provided default
        # to 10**nlog<0.05
        sig = sigs.get(
            k,
            (10 ** (-curr_nlog_dsea_res)) < 0.05
        )
        if any(sig):
            sig_nlog_dsea_res = curr_nlog_dsea_res.copy()
            sig_nlog_dsea_res[~sig] = np.NaN
            ax.scatter(
                x=np.arange(len(sig_nlog_dsea_res)) + offsets[k],
                y=sig_nlog_dsea_res + 1,
                s=10,
                marker=r"$*$",
                c=".3",
            )
    ax.set_ylabel("-log(p)")
    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1] + 1)
    ax.set_xlabel("")
    ax.set_xticklabels([s.replace(" ", "\n") for s in DSEA_STAGES], fontsize=9)
    ax.legend(**legend_kwargs)
    return ax

def dsea_spin(
    X, parcellation_name, memmap, n_perm=1000,
    seed=0, mirror='new', n_genes=500, pSI='0.05',
    fdr=True
):
    """
    Performs spin permutation test on dSEA

    Parameters
    ---------
    X: (pd.DataFrame)
    parcellation_name: (str)
    memmap: (str)
        path to memmap file to store null distribution
    n_perm: (int)
        Number of spin permutations
    seed: (int)
        PLS random seed
    mirror: {'new', 'old', str}
        CSEA tool url
    n_genes: (int)
        Number of top genes to select
    pSI: (str)
        specificity index threshold
    fdr: (bool)
        apply FDR correction on p-values

    Returns
    -------
    test_res: (pd.Series)
        test results
    p_vals: (pd.Series)
        spin p-values
    raw_dsea_res: (dict)
        raw dSEA results
    top_genes: (list of dict)
        top genes from PLS with 'pos' and 'neg' weights
    pls: (sklearn.cross_decomposition.PLSRegression)
        fitted PLS model
    """
    # run dSEA on the true map
    test_res, raw_dsea_res, top_genes, pls = dsea(
        X, parcellation_name, 
        seed=seed, n_genes=n_genes,
        mirror=mirror, pSI=pSI, fdr=fdr
    )
    # create spin surrogates
    rotated_parcels = get_rotated_parcels(
        parcellation_name, n_perm, return_indices=False, seed=seed,
    )
    # create (or continue) null distribution as a memmap
    # (given each iteration takes a long time and may fail
    # due to connection issues, this will prevent waste of
    # time in case an error occurs during iterations)
    if os.path.exists(memmap):
        null_dist = np.memmap(memmap, dtype='float64', mode='r+', shape=(test_res.shape[0], n_perm))
        # determine starting permutation as number
        # of permutations in which first row is not NaN
        start_perm = np.sum(~np.isnan(null_dist[0]))
        print(f"continuing from permutation {start_perm}")
    else:
        null_dist = np.memmap(memmap, dtype='float64', mode='w+', shape=(test_res.shape[0], n_perm))
        # set all to NaNs which is needed to determine
        # start_perm in the next retries
        null_dist[:, :] = np.NaN
        start_perm = 0 
    for i in tqdm(range(start_perm, n_perm)):
        # create surrogate data with rotated parcels
        surrogate_data = X.loc[rotated_parcels[:, i]]
        surrogate_data.index = X.index
        # run dSEA on the surrogate data
        null_res, _, _, _ = dsea(
            surrogate_data, parcellation_name, 
            seed=seed, n_genes=n_genes,
            mirror=mirror, pSI=pSI, fdr=fdr
        )
        null_dist[:, i] = null_res.values
    # p-value
    p_vals = (np.abs(null_dist) > np.abs(test_res.values[:, np.newaxis])).mean(axis=1)
    p_vals = pd.Series(p_vals, index=test_res.index)
    return test_res, p_vals, raw_dsea_res, top_genes, pls