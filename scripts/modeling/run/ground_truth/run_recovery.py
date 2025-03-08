import time
import numpy as np
from cubnm import optimize, _setup_opts
import os
import json
import sys

RESULTS_DIR = os.path.dirname(os.path.join(os.environ["PROJECT_DIR"], "results"))
INPUT_DIR = os.environ["INPUT_DIR"]
OUTPUT_DIR = os.environ["OUTPUT_DIR"]

# ensure cubnm is installed with the CUBNM_NOISE_WHOLE env variable set to 1
assert not _setup_opts.noise_segment_flag, (
    "cubnm must be installed with the CUBNM_NOISE_WHOLE env variable set to 1"
)

def run_recovery(sim_seed, n_runs=2):
    """
    Runs the recovery CMA-ES optimiation

    Parameters
    ----------
    sim_seed : int
        Random seed for the simulation
    n_runs : int
        Number of optimization runs
    """
    maps_path = os.path.join(INPUT_DIR, "6maps_schaefer-100_zscore.txt")
    sc_path = os.path.join(
        INPUT_DIR,
        "micamics",
        "SC",
        "group-all",
        "ctx_parc-schaefer-100_approach-median_mean001_desc-strength.txt",
    )
    fc_path = os.path.join(RESULTS_DIR, "ground_truth", "synthetic", "fc_tril.txt")
    fcd_path = os.path.join(RESULTS_DIR, "ground_truth", "synthetic", "fcd_tril.txt")
    opt_seeds = list(range(1, n_runs + 1))
    problems = []
    optimizers = []
    for opt_seed in opt_seeds:
        out_path = os.path.join(RESULTS_DIR, "ground_truth", "recovery", f"sim-{sim_seed}_opt-{opt_seed}")
        if os.path.exists(out_path):
            if len(os.listdir(out_path)) > 0:
                print(f"Output directory {out_path} already exists and is not empty. Skipping.")
                continue
        # set up problem objects for all cases
        problem = optimize.BNMProblem(
            model='rWW',
            params={
                "G": (0.5, 4.0),
                "wEE": (0.05, 0.75),
                "wEI": (0.05, 0.75),
            },
            het_params=['wEE', 'wEI'],
            max_fic_trials=10,
            duration=450,
            TR=3,
            window_size=10,
            window_step=2,
            sc=sc_path,
            emp_fc_tril=np.loadtxt(fc_path),
            emp_fcd_tril=np.loadtxt(fcd_path),
            gof_terms=['+fc_corr', '-fc_diff', '-fcd_ks'],
            maps=maps_path,
            out_dir=out_path,
            bw_params="heinzle2016-3T",
            exc_interhemispheric=True,
            rand_seed=sim_seed,
        )
        problems.append(problem)
        optimizer = optimize.CMAESOptimizer(
            popsize=210,
            n_iter=81,
            seed=opt_seed,
            algorithm_kws=dict(tolfun=5e-3),
        )
        optimizers.append(optimizer)
    start = time.time()
    optimizers = optimize.batch_optimize(optimizers, problems)
    print(
        f"CMAES for {len(problems)} runs with 210 particles and 81 iterations took a total walltime of {time.time()-start}s"
    )


if __name__ == "__main__":
    sim_seed = int(sys.argv[1])
    run_recovery(sim_seed)