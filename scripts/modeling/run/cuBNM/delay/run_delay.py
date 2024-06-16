"""
Given a CMAES result with a given simulation seed, runs the simulation using the best
parameters with alternative simulation seeds
It is important to install cuBNM with the CUBNM_NOISE_WHOLE env variable set to 1
"""

import os
import sys
import numpy as np
import gc

from cuBNM import sim

TRs = {"pnc": 3, "imagen": 2.2}


def run_delay(dataset, cmaes_log_path, sc_path, fc_tril_path):
    """
    Gets the best parameters from a CMAES log file and runs the simulation with
    delay

    Parameters
    ----------
    cmaes_log_path : str
        Path to the CMAES log file created by `bnm` repository
        using the typical config

    """
    # use the following velocities
    v_list = [1, 2, 3, 4, 5, 6]
    if os.path.exists(
        os.path.join(
            cmaes_log_path.replace(".txt", "_cuBNM"), f"v-{v_list[-1]}", "score.csv"
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
    # get wEE and wEI arrays from simulation output
    sims_out_prefix = cmaes_log_path.replace("/cmaes", "/sims").replace(".txt", "")
    param_lists["wEE"] = np.loadtxt(sims_out_prefix + "_w_EE.txt")[np.newaxis, :]
    param_lists["wEI"] = np.loadtxt(sims_out_prefix + "_w_EI.txt")[np.newaxis, :]

    # load empirical FC and FCD
    emp_fc_tril = np.loadtxt(fc_tril_path)
    emp_fcd_tril = np.loadtxt(fc_tril_path.replace("FCtril", "FCDtril"))

    # run original simulation
    # with and without sync_msec set to True
    for sync_msec in [False, True]:
        print(f"Running original simulation with no delay (sync_msec = {sync_msec})")
        orig_sg = sim.rWWSimGroup(
            duration=450,
            TR=TRs[dataset],
            sc_path=sc_path,
            sc_dist_path=None,
            bw_params="heinzle2016-3T",
            out_dir=os.path.join(
                cmaes_log_path.replace(".txt", "_cuBNM"),
                f"v-orig_syncmsec-{int(sync_msec)}",
            ),
            max_fic_trials=10,
            sim_verbose=True,
            progress_interval=1000,
        )
        orig_sg.N = 1
        orig_sg.param_lists["G"] = param_lists["G"]
        orig_sg.param_lists["wEE"] = param_lists["wEE"]
        orig_sg.param_lists["wEI"] = param_lists["wEI"]
        orig_sg.sync_msec = sync_msec
        orig_sg.run()
        orig_sg.save()
        orig_sg.score(emp_fc_tril, emp_fcd_tril).to_csv(
            os.path.join(orig_sg.out_dir, "score.csv")
        )
        del orig_sg
        gc.collect()

    # run and score simulations with alternative velocities
    for v in v_list:
        print("Running velocity", v)
        sg = sim.rWWSimGroup(
            duration=450,
            TR=TRs[dataset],
            sc_path=sc_path,
            sc_dist_path=sc_path.replace("strength", "length"),
            bw_params="heinzle2016-3T",
            out_dir=os.path.join(cmaes_log_path.replace(".txt", "_cuBNM"), f"v-{v}"),
            max_fic_trials=10,
            sim_verbose=True,
            progress_interval=1000,
        )
        sg.N = 1
        sg.param_lists["G"] = param_lists["G"]
        sg.param_lists["wEE"] = param_lists["wEE"]
        sg.param_lists["wEI"] = param_lists["wEI"]
        sg.param_lists["v"] = np.array([float(v)])
        sg.run(force_reinit=True)
        sg.save()
        sg.score(emp_fc_tril, emp_fcd_tril).to_csv(
            os.path.join(sg.out_dir, "score.csv")
        )
        del sg
        gc.collect()


if __name__ == "__main__":
    dataset = sys.argv[1]
    cmaes_log_path = sys.argv[2]
    sc_path = sys.argv[3]
    fc_tril_path = sys.argv[4]
    run_delay(dataset, cmaes_log_path, sc_path, fc_tril_path)
