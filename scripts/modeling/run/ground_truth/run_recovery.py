import time
import numpy as np
from cubnm import optimize, _setup_opts
import os
import json
import sys

EXP_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.environ["INPUT_DIR"]
OUTPUT_DIR = os.environ["OUTPUT_DIR"]

# ensure cubnm is installed with the CUBNM_NOISE_WHOLE env variable set to 1
assert not _setup_opts.noise_segment_flag, (
    "cubnm must be installed with the CUBNM_NOISE_WHOLE env variable set to 1"
)

def run_recovery(sim_seed, opt_seed):
    """
    Runs the recovery CMA-ES optimiation
    """
    maps_path = os.path.join(INPUT_DIR, "6maps_schaefer-100_zscore.txt")
    sc_path = os.path.join(
        INPUT_DIR,
        "micamics",
        "SC",
        "group-all",
        "ctx_parc-schaefer-100_approach-median_mean001_desc-strength.txt",
    )
    fc_path = os.path.join(EXP_DIR, "synthetic", "fc_tril.txt")
    fcd_path = os.path.join(EXP_DIR, "synthetic", "fcd_tril.txt")
    out_path = os.path.join(EXP_DIR, "recovery", f"sim-{sim_seed}_opt-{opt_seed}")
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
    optimizer = optimize.CMAESOptimizer(
        popsize=210,
        n_iter=81,
        seed=opt_seed,
        algorithm_kws=dict(tolfun=5e-3),
    )
    optimizer.setup_problem(problem)
    optimizer.optimize()
    optimizer.save()



if __name__ == "__main__":
    sim_seed = int(sys.argv[1])
    opt_seed = int(sys.argv[2])
    run_recovery(sim_seed, opt_seed)