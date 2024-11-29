import sys
import os
import copy
import numpy as np
import pandas as pd
import scipy.stats

PROJECT_DIR = os.environ.get("PROJECT_DIR")
CODE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.append(CODE_DIR)
from proc_rest import post_fmriprep


def check_fit_bold(
    emp_fc_tril,
    emp_fcd_tril,
    sims_dir,
    sim_bold_filename,
    dataset,
    sim_vol_remove=10,
    return_fc_fcd=False,
    exc_interhemispheric=True,
):
    """
    Checks fit of empirical and simulated FC. Simulations are
    read from `sims_dir` which includes files generated using
    https://github.com/amnsbr/bnm_cuda


    Parameters
    ---------
    emp_fc_tril: (np.ndarray)
        lower triangle of empirical FC
    emp_fcd_tril: (np.ndarray)
        lower triangle of empirical FCD
    sims_dir: (pathlike str)
        location of the simulation files
    sim_bold_filename: (str)
        filename of the BOLD text file
    dataset: (str)
    sim_vol_remove: (int)
        number of BOLD volumes to remove from the start of simulation
    return_fc_fcd: (bool)
        returns simulated fc and fcd in addition to fit results
    exc_interhemispheric: (bool)
        exclude interhemispheric connections from simulated FC and FCD

    Returns
    ------
    res: (dict)
        with keys 'G', 'wee', 'wei', 'wie_scale', 'corr', 'ks', 'cost'
    if `return_fc_fcd`, will additionally return:
    sim_fc_tril, sim_fcd_tril: (np.ndarray)
    """
    res = {}
    # load sim and remove volumes
    # using pd.read_csv to handle occasional kernel crashing
    # with np.loadtxt
    sim_bold = pd.read_csv(
            os.path.join(sims_dir, sim_bold_filename),
            delim_whitespace=True,
            header=None
        ).squeeze().values
    if sim_bold.size == 0:
        print(sim_bold_filename, "is empty")
        return res
    sim_bold = sim_bold[sim_vol_remove:, :]
    # FC correlation and KS
    sim_fc = np.corrcoef(sim_bold.T)
    rh_idx = int(
        sim_fc.shape[0] / 2
    )  # assumes symmetric number of parcels ordered L->R
    if exc_interhemispheric:
        sim_fc[:rh_idx, rh_idx:] = np.NaN
        sim_fc[rh_idx:, :rh_idx] = np.NaN
    sim_fc_tril = sim_fc[np.tril_indices_from(sim_fc, -1)]
    sim_fc_tril = sim_fc_tril[~np.isnan(sim_fc_tril)]
    res["fc_corr"] = np.corrcoef(emp_fc_tril, sim_fc_tril)[0, 1]
    res["fc_diff"] = np.abs(np.nanmean(sim_fc_tril) - np.nanmean(emp_fc_tril))
    # FCD KS
    sim_bold_df = pd.DataFrame(sim_bold.T)
    hemi_parcs = {"L": list(range(rh_idx)), "R": list(range(rh_idx, sim_fc.shape[0]))}
    sim_fcd, sim_fc_windows = post_fmriprep.calculate_fcd(
        sim_bold_df,
        post_fmriprep.WINDOW_SIZE[dataset],
        post_fmriprep.STEP[dataset],
        exc_interhemispheric=exc_interhemispheric,
        hemi_parcs=hemi_parcs,
    )
    sim_fcd_tril = sim_fcd[np.tril_indices_from(sim_fcd, -1)]
    res["fcd_ks"] = scipy.stats.ks_2samp(emp_fcd_tril, sim_fcd_tril).statistic
    if return_fc_fcd:
        return res, sim_fc_tril, sim_fcd_tril
    else:
        return res


def load_cmaes(
    cmaes_dir,
    sims_dir,
    params,
    het_params,
    emp_fc_tril,
    emp_fcd_tril,
    n_vols_remove,
    dataset,
    sc_path=None,
    itMax=81,
    lmbda=210,
    SeedMW=1,
    SeedSim=410,
    fc_prefix="ctx_parc-schaefer-100_hemi-LR_highpass-013_lowpass-none_exc-inter",
    nodes=100,
    exc_interhemispheric=True,
):
    """
    Loads the goodness-of-fit and regional variables of the optimal
    point of a CMAES run done using https://github.com/amnsbr/bnm_cuda

    Parameters
    ----------
    cmaes_dir: (str)
    sims_dir: (str)
    params: (dict)
        with keys "G", "wee", "wei", and "wie"
    het_params: (str)
        heterogeneity parameters separated by "-"
    emp_fc_tril: (np.ndarray)
    emp_fcd_tril: (np.ndarray)
    n_vols_remove: (int)
        number of BOLD volumes to remove from the start of simulation
    dataset: (str)
    sc_path: (str)
        if included, sc_fc coupling will be calculated
    itMax: (int)
        number of CMAES iterations
    lmbda: (int)
        number of CMAES particles
    SeedMW: (int)
        random seed of CMAES
    SeedSim: (int)
        random seed of simulation
    fc_prefix: (str)
        FC(D) filename prefix
    nodes: (int)
        nubmer of nodes
    exc_interhemispheric: (bool)
        exclude interhemispheric connections from simulated FC and FCD

    Returns
    -------
    res: (dict of float/int)
        goodness-of-fit and related measures of the optima
    regional_vars: (dict of np.ndarray)
        regional variables including regional parameters and
        mean state variables
    """
    params = copy.deepcopy(params)
    prefix = f'{fc_prefix}_G_{params["G"]}_wee_{params["wee"]}_wei_{params["wei"]}_wie_{params["wie"]}_het-{het_params}_SeedMW-{SeedMW}_SeedSim-{SeedSim}_n-{itMax}x{lmbda}'
    cmaes_log_path = os.path.join(cmaes_dir, f"{prefix}.txt")
    if not (os.path.exists(cmaes_log_path)):
        print(f"{cmaes_log_path} does not exist")
        return
    best_params_str = ""
    with open(cmaes_log_path, "r") as f:
        lines = f.read().split("\n")
        for line in lines:
            if line.startswith("Iteration"):
                last_iteration = int(line.split(" ")[1])
            elif line.startswith("Best"):
                best_params_str = line
    if not best_params_str:
        print(f"{cmaes_log_path} is incomplete")
        return
    best_params_str = best_params_str.replace("Best feasible parameters: ", "")
    best_params_str = best_params_str.rstrip(",")
    best_params = {
        s.split("=")[0]: float(s.split("=")[1]) for s in best_params_str.split(",")
    }
    # GOF calculations
    best_params["fic_penalty"] = float(best_params["fic_penalty"])
    best_params["gof_fic_penalty"] = float(best_params["gof_fic_penalty"])
    best_params["gof_printed"] = float(best_params["gof"])
    res = copy.deepcopy(best_params)
    res["last_iteration"] = last_iteration
    if sc_path:
        sc = np.loadtxt(sc_path)
        if exc_interhemispheric:
            rh_idx = int(sc.shape[0] / 2)
            sc[:rh_idx, rh_idx:] = np.NaN
            sc[rh_idx:, :rh_idx] = np.NaN
        sc_tril = sc[np.tril_indices_from(sc, -1)]
        sc_tril = sc_tril[~np.isnan(sc_tril)]
        res["sc_fc"] = np.corrcoef(sc_tril, emp_fc_tril)[0, 1]
    best_sim_file = prefix + "_bold.txt"
    fit_py = check_fit_bold(
        emp_fc_tril,
        emp_fcd_tril,
        sims_dir,
        best_sim_file,
        dataset,
        sim_vol_remove=n_vols_remove,
        exc_interhemispheric=exc_interhemispheric,
    )
    res["fc_corr"] = fit_py["fc_corr"]
    res["fcd_ks"] = fit_py["fcd_ks"]
    res["fc_diff"] = fit_py["fc_diff"]
    res["cost"] = -res["fc_corr"] + res["fc_diff"] + res["fcd_ks"]
    res["gof"] = -res["cost"]
    res["fic_penalty"] = res["gof"] - res["gof_fic_penalty"]
    # fic adjustment
    res["fic_failed"] = 1 in np.loadtxt(
        os.path.join(sims_dir, f"{prefix}_fic_failed.txt")
    )
    with open(os.path.join(sims_dir, f"{prefix}_fic_ntrials.txt"), "r") as f:
        res["fic_ntrials"] = int(f.read().split(" ")[0])
    # parameters and average states
    regional_vars = {}
    for var in ["I_E", "r_E", "S_E", "I_I", "r_I", "S_I", "w_EE", "w_EI", "w_IE"]:
        regional_vars[var] = pd.read_csv(
                os.path.join(sims_dir, f"{prefix}_{var}.txt"),
                delim_whitespace=True,
                header=None,
            ).squeeze().values
    regional_vars = pd.DataFrame(regional_vars)
    return res, regional_vars

def load_cmaes_hist(
    cmaes_dir,
    params,
    het_params,
    itMax=81,
    lmbda=210,
    SeedMW=1,
    fc_prefix="ctx_parc-schaefer-100_hemi-LR_highpass-013_lowpass-none_exc-inter",
):
    """
    Loads the history of particles (simulations) done during
    a CMAES run done using https://github.com/amnsbr/bnm_cuda

    Parameters
    ----------
    cmaes_dir: (str)
    params: (dict)
        with keys "G", "wee", "wei", and "wie"
    het_params: (str)
        heterogeneity parameters separated by "-"
    itMax: (int)
        number of CMAES iterations
    lmbda: (int)
        number of CMAES particles
    SeedMW: (int)
        random seed of CMAES
    fc_prefix: (str)
        FC(D) filename prefix

    Returns
    -------
    cmaes: (pd.DataFrame)
        parameters and goodness-of-fit measures of all particles
    optima: (pd.DataFrame)
        parameters and goodness-of-fit measures of optima per generation
    optima_cumulative: (pd.DataFrame)
        parameters and goodness-of-fit measures of optima up to each generation
    """
    params = copy.deepcopy(params)
    prefix = f'{fc_prefix}_G_{params["G"]}_wee_{params["wee"]}_wei_{params["wei"]}_wie_{params["wie"]}_het-{het_params}_SeedMW-{SeedMW}_SeedSim-410_n-{itMax}x{lmbda}'
    cmaes_log_path = os.path.join(cmaes_dir, f"{prefix}.txt")
    if not (os.path.exists(cmaes_log_path)):
        print(f"{cmaes_log_path} does not exist")
        return
    with open(cmaes_log_path, "r") as f:
        lines = f.read().split("\n")
    cmaes = []
    for line in lines:
        if line.startswith("Iteration"):
            curr_iteration = int(line.split(" ")[1])
        elif line.startswith("\t"):
            curr_sample = dict(
                [p.strip().split("=") for p in line[line.find("G") :].split(",")]
            )
            curr_sample["iteration"] = curr_iteration
            cmaes.append(curr_sample)
    cmaes = pd.DataFrame(cmaes).astype("float")
    optima_idx = cmaes.groupby("iteration")["gof"].idxmax().astype("int")
    optima = cmaes.loc[optima_idx.values]
    optima_cumulative = pd.Series(index=optima.index)
    for it in optima_cumulative.index:
        optima_cumulative.loc[it] = cmaes.loc[cmaes["iteration"] <= it, "gof"].max()
    return cmaes, optima, optima_cumulative
