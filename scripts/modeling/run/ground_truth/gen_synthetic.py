import os
import numpy as np

from cubnm import sim, _setup_opts

RESULTS_DIR = os.path.dirname(os.path.join(os.environ["PROJECT_DIR"], "results"))
INPUT_DIR = os.environ["INPUT_DIR"]

# ensure cubnm is installed with the CUBNM_NOISE_WHOLE env variable set to 1
assert not _setup_opts.noise_segment_flag, (
    "cubnm must be installed with the CUBNM_NOISE_WHOLE env variable set to 1"
)

def run():
    """
    Runs a simulation which will be used as a synthetic ground truth
    to test parameter recovery
    """
    # use micamics group-averaged SC
    sc_path = os.path.join(
        INPUT_DIR, 'micamics', 'SC', 'group-all',
        'ctx_parc-schaefer-100_approach-median_mean001_desc-strength.txt'
    )
    # set parameters (roughly approximate optimal parameters
    # of the second age group to ensure parameters are in
    # a reasonable range)
    param_lists = {}
    param_lists["G"] = np.array([2.5])
    # set wEE and wEI based on maps and known base and coefficients
    # this code is similar to the code used in 
    # `cuBNM.optimize.BNMProblem._set_sim_params`
    maps = np.loadtxt(os.path.join(INPUT_DIR, '6maps_schaefer-100_zscore.txt'))
    het_params = ['wEE', 'wEI']
    biases = {
        'wEE': 0.075,
        'wEI': 0.40,
    }
    param_lists["wEE"] = np.tile(biases['wEE'], (1, 100))
    param_lists["wEI"] = np.tile(biases['wEI'], (1, 100))
    coefficients = {
        'wEEscale0': 0.15,
        'wEIscale0': 0.05,
        'wEEscale1': -0.30,
        'wEIscale1': 0.15,
        'wEEscale2': -0.05,
        'wEIscale2': -0.15,
        'wEEscale3': 0.15,
        'wEIscale3': 0.05,
        'wEEscale4': 0.01,
        'wEIscale4': 0.35,
        'wEEscale5': 0.20,
        'wEIscale5': 0.20,
    }
    for param in het_params:
        param_scalers = np.ones(100)
        for map_idx in range(maps.shape[0]):
            param_scalers += (
                maps[map_idx, :] * coefficients[f"{param}scale{map_idx}"]
            )
        param_lists[param] *= param_scalers[np.newaxis, :]
        # shift if there are values < 0.001 (incl. negative values)
        if param_lists[param][0, :].min() < 0.001:
            param_lists[param][0, :] -= (
                param_lists[param][0, :].min()
                - 0.001
            )
    # run the simulation
    sg = sim.rWWSimGroup(
        duration=450,
        TR=3,
        window_size=10,
        window_step=2,
        sc=sc_path,
        bw_params="heinzle2016-3T",
        exc_interhemispheric=True,
        out_dir=os.path.join(RESULTS_DIR, "ground_truth", "synthetic"),
        max_fic_trials=10,
        rand_seed=410,
    )
    sg.N = 1
    for param in ["G", "wEE", "wEI"]:
        sg.param_lists[param] = param_lists[param].copy()
    sg.run()
    sg.save()
    np.savetxt(
        os.path.join(sg.out_dir, "fc_tril.txt"),
        sg.sim_fc_trils[0]
    )
    np.savetxt(
        os.path.join(sg.out_dir, "fcd_tril.txt"),
        sg.sim_fcd_trils[0]
    )


if __name__ == "__main__":
    run()
