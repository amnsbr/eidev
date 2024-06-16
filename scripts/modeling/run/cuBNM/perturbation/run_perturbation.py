"""
Given a CMAES result with a given simulation seed, runs the simulation using the best
parameters with perturbed parameters by -10% or 10%
It is important to install cuBNM with the CUBNM_NOISE_WHOLE env variable set to 1
"""

import os
import sys
import numpy as np

from cuBNM import sim

TRs = {"pnc": 3, "imagen": 2.2}


def run(dataset, cmaes_log_path, sc_path, fc_tril_path):
    """
    Gets the best parameters from a CMAES log file and runs the simulation with
    altered parameters in range of -10% to 10%

    Parameters
    ----------
    cmaes_log_path : str
        Path to the CMAES log file created by `bnm` repository
        using the typical config
    """
    if os.path.exists(
        os.path.join(
            cmaes_log_path.replace(".txt", "_cuBNM"),
            "int-G_x1.10",
            "score.csv",
        )
    ):
        print("Already done. Exiting...")
        return
    # get parameters from CMAES log and associated simulation outputs
    param_lists = {}
    # get G from CMAES log best parameters
    with open(cmaes_log_path, "r") as f:
        lines = f.read().split("\n")
    for line in lines:
        if line.startswith("Best"):
            best_params_str = line
    if not best_params_str:
        print(f"{cmaes_log_path} is incomplete")
        return
    best_params_str = best_params_str.replace("Best feasible parameters: ", "")
    best_params = {
        s.split("=")[0]: float(s.split("=")[1]) for s in best_params_str.split(",")
    }
    param_lists["G"] = np.array([best_params["G"]])
    # get wEE, wEI and FIC wIE arrays from simulation output
    sims_out_prefix = cmaes_log_path.replace("/cmaes", "/sims").replace(".txt", "")
    param_lists["wEE"] = np.loadtxt(sims_out_prefix + "_w_EE.txt")[np.newaxis, :]
    param_lists["wEI"] = np.loadtxt(sims_out_prefix + "_w_EI.txt")[np.newaxis, :]
    param_lists['wIE'] = np.loadtxt(sims_out_prefix + '_w_IE.txt')[np.newaxis,:]

    # run the altered simulations (without FIC)
    for alt_param in ["wEE", "wEI", "wIE", "G"]:
        for ratio in [0.9, 1.1]:
            print(f"Running simulation with {alt_param}x{ratio}")
            sg = sim.rWWSimGroup(
                duration=450,
                TR=TRs[dataset],
                sc_path=sc_path,
                bw_params="heinzle2016-3T",
                out_dir=os.path.join(
                    cmaes_log_path.replace(".txt", "_cuBNM"), f"int-{alt_param}_x{round(ratio, 2)}"
                ),
                do_fic=False,
                max_fic_trials=0,
                extended_output=True,
                extended_output_ts=True,
                rand_seed=410,
            )
            sg.N = 1
            for curr_param in ["wEE", "wEI", "wIE", "G"]:
                sg.param_lists[curr_param] = param_lists[curr_param].copy()
                if curr_param == alt_param:
                    sg.param_lists[curr_param] *= ratio
            sg.fic_unstable = False # workaround for a bug in .save()
            sg.fic_failed = False # "
            sg.fic_ntrials = 0 # "
            sg.run()
            sg.save()


if __name__ == "__main__":
    dataset = sys.argv[1]
    cmaes_log_path = sys.argv[2]
    sc_path = sys.argv[3]
    fc_tril_path = sys.argv[4]
    run(dataset, cmaes_log_path, sc_path, fc_tril_path)
