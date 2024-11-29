import os
import numpy as np
import pandas as pd
from cubnm import sim, _setup_opts

PROJECT_DIR = os.environ['PROJECT_DIR']
INPUT_DIR = os.environ['INPUT_DIR']
PNC_OUTPUT_DIR = os.path.join(os.environ["PNC_PROJECT_DIR"], "output")
RESULTS_DIR = os.path.join(PROJECT_DIR, 'results')

# ensure cubnm is installed with the CUBNM_NOISE_WHOLE env variable set to 1
assert not _setup_opts.noise_segment_flag, (
    "cubnm must be installed with the CUBNM_NOISE_WHOLE env variable set to 1"
)

def run_perturbation_all(orig_seed=410, n_seeds=50):
    """
    Given the optimal simulations of 40 randomly selected subjects,
    performs perturbed simulations with each parameter altered by +10% or -10%
    """
    # load optimal simulations data of all subjects
    optima_csv = os.path.join(RESULTS_DIR, 'pnc_optima.csv')
    if not os.path.exists(optima_csv):
        print("Run Figure 2 first")
        return
    optima = pd.read_csv(optima_csv, index_col=0)
    # select best runs
    best_runs = optima.groupby('sub')['gof'].idxmax().values
    data = optima.loc[best_runs].set_index('sub')
    # randomly select 40 subjects
    data = data.sample(n=40, random_state=0)
    print(data.index)
    # curate optimal (original) parameters and SCs
    param_lists = {
        'G':[],
        'wEE': [],
        'wEI': [],
        'wIE': [],
    }
    SCs = []    
    for sub, row in data.iterrows():
        SeedMW = row['SeedMW']
        cmaes_log_path = os.path.join(
            PNC_OUTPUT_DIR, 'sim', sub, 'ctx_parc-schaefer-100_mean001_thresh-1',
            '6maps_schaefer-100_zscore', 'cmaes_multimaps_gpu', 
            f'ctx_parc-schaefer-100_hemi-LR_highpass-013_lowpass-none_exc-inter_G_0.5-4_wee_0.05-0.75_wei_0.05-0.75_wie_0_het-wee-wei_SeedMW-{SeedMW}_SeedSim-{orig_seed}_n-81x210.txt'
        )
        sc_path = os.path.join(
            PNC_OUTPUT_DIR, 'SC', sub, 
            'ctx_parc-schaefer-100_mean001_thresh-1_desc-strength.txt'
        )    
        SCs.append(np.loadtxt(sc_path))
        # get parameters from CMAES log and associated simulation outputs
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
        param_lists["G"].append(best_params["G"])
        # get wEE, wEI and wIE arrays from simulation output
        sims_out_prefix = cmaes_log_path.replace("/cmaes", "/sims").replace(".txt", "")
        param_lists["wEE"].append(np.loadtxt(sims_out_prefix + "_w_EE.txt"))
        param_lists["wEI"].append(np.loadtxt(sims_out_prefix + "_w_EI.txt"))
        param_lists["wIE"].append(np.loadtxt(sims_out_prefix + "_w_IE.txt"))
    # concatenate parameters of all subjects
    param_lists = {k:np.array(v) for k,v in param_lists.items()}
    # concatenate SCs
    SCs = np.array(SCs)
    # define sim group
    # skip doing FIC since wIE is set to the values based on FIC
    # solution of optimal simulations, and we'd want to isolate the
    # effect of perturbations in each parameter without changing
    # any other parameters (which would happen if FIC is on)
    sg = sim.rWWSimGroup(
        duration=450,
        TR=3,
        window_size=10,
        window_step=2,
        sc=SCs[0],
        bw_params="heinzle2016-3T",
        exc_interhemispheric=True,
        gof_terms=['+fc_corr', '-fc_diff', '-fcd_ks'],
        do_fic=False,
        max_fic_trials=0,
        states_ts=True,
        sim_verbose=True,
        rand_seed=410,
    )
    sg.N = len(param_lists['G'])
    # the current version of cubnm does not support
    # variable SCs assigned in SimGroup initialization
    # and they must be assigned after it's instantiated
    sg.sc = SCs
    sg.sc_indices = np.arange(sg.N)
    sg.param_lists["G"] = param_lists["G"]
    sg.param_lists["wEE"] = param_lists["wEE"]
    sg.param_lists["wEI"] = param_lists["wEI"]
    sg.param_lists["wIE"] = param_lists["wIE"]
    # run non-altered simulations
    sg.out_dir = os.path.join(
        RESULTS_DIR, "perturbation",
        "orig",
    )
    os.makedirs(sg.out_dir, exist_ok=True)
    sg.run()
    sg.save()
    # run the perturbed simulations
    for alt_param in ["wEE", "wEI", "wIE", "G"]:
        for ratio in [0.9, 1.1]:
            print(f"Running simulation with {alt_param}x{ratio}")
            sg.out_dir = os.path.join(
                RESULTS_DIR, "perturbation",
                f"{alt_param}_x{round(ratio, 2)}",
            )
            if os.path.exists(os.path.join(sg.out_dir, 'scores.csv')):
                print("already done")
                continue
            os.makedirs(sg.out_dir, exist_ok=True)
            # alter a parameter while keeping the rest fixed to optimal
            for curr_param in ["wEE", "wEI", "wIE", "G"]:
                sg.param_lists[curr_param] = param_lists[curr_param].copy()
                if curr_param == alt_param:
                    sg.param_lists[curr_param] *= ratio
            # run and save
            sg.run()
            sg.it = 1 # to have it1.npz as filenames
            sg.save()

if __name__ == "__main__":
    run_perturbation_all()
