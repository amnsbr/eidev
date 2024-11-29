import os
import numpy as np
import pandas as pd
import scipy
from tqdm import tqdm
from cubnm import sim, _setup_opts

PROJECT_DIR = os.environ['PROJECT_DIR']
INPUT_DIR = os.environ['INPUT_DIR']
PNC_OUTPUT_DIR = os.path.join(os.environ["PNC_PROJECT_DIR"], "output")
RESULTS_DIR = os.path.join(PROJECT_DIR, 'results')

# ensure cubnm is installed with the CUBNM_NOISE_WHOLE env variable set to 1
assert not _setup_opts.noise_segment_flag, (
    "cubnm must be installed with the CUBNM_NOISE_WHOLE env variable set to 1"
)

def run_sim_seed_all(orig_seed=410, n_seeds=50):
    """
    Gets the best parameters from a CMAES log file and runs the simulation with
    alternative seeds
    """
    # load optimal simulations data of all subjects
    optima_csv = os.path.join(PROJECT_DIR, 'results', 'pnc_optima.csv')
    if not os.path.exists(optima_csv):
        print("Run Figure 2 first")
        return
    optima = pd.read_csv(optima_csv, index_col=0)
    # limit the optima to subsample
    subsample = np.loadtxt(os.path.join(INPUT_DIR, 'pnc_subsample_200.txt'), dtype=str)
    optima = optima.loc[optima['sub'].isin(subsample)]
    # select best runs
    best_runs = optima.groupby('sub')['gof'].idxmax().values
    data = optima.loc[best_runs].set_index('sub')
    param_lists = {
        'G':[],
        'wEE': [],
        'wEI': [],
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
        # get wEE and wEI arrays from simulation output
        sims_out_prefix = cmaes_log_path.replace("/cmaes", "/sims").replace(".txt", "")
        param_lists["wEE"].append(np.loadtxt(sims_out_prefix + "_w_EE.txt"))
        param_lists["wEI"].append(np.loadtxt(sims_out_prefix + "_w_EI.txt"))
    # concatenate parameters of all subjects
    param_lists = {k:np.array(v) for k,v in param_lists.items()}
    # concatenate SCs
    SCs = np.array(SCs)
    # define sim group
    sg = sim.rWWSimGroup(
        duration=450,
        TR=3,
        window_size=10,
        window_step=2,
        sc=SCs[0],
        bw_params="heinzle2016-3T",
        exc_interhemispheric=True,
        gof_terms=['+fc_corr', '-fc_diff', '-fcd_ks'],
        max_fic_trials=10,
        sim_verbose=True,
        rand_seed=orig_seed,
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
    for rand_seed in [orig_seed]+list(range(n_seeds)):
        print("Running original simulation with seed", rand_seed)
        sg.rand_seed = rand_seed
        sg.out_dir = os.path.join(
            RESULTS_DIR, "sim_seed",
            f"SeedSim-{int(rand_seed)}",
        )
        if os.path.exists(os.path.join(sg.out_dir, 'scores.csv')):
            print("already done")
            continue
        os.makedirs(sg.out_dir, exist_ok=True)
        sg.run()
        # sg.save() # this ends up being very large!
        # so just save what we need
        np.savez_compressed(os.path.join(sg.out_dir, 'states.npz'), **sg.sim_states)
        # check fit
        scores = {}
        for i, sub in tqdm(enumerate(data.index)):
            fc_tril_path = os.path.join(
                PNC_OUTPUT_DIR, 'FC', sub,
                'ctx_parc-schaefer-100_hemi-LR_highpass-013_lowpass-none_exc-inter_desc-FCtril.txt'
            )
            # load empirical FC and FCD
            emp_fc_tril = np.loadtxt(fc_tril_path)
            emp_fcd_tril = np.loadtxt(fc_tril_path.replace("FCtril", "FCDtril"))
            # score
            scores[sub] = {
                '+fc_corr': scipy.stats.pearsonr(sg.sim_fc_trils[i], emp_fc_tril).statistic,
                '-fc_diff': -np.abs(sg.sim_fc_trils[i].mean() - emp_fc_tril.mean()),
                '-fcd_ks': -scipy.stats.ks_2samp(sg.sim_fcd_trils[i], emp_fcd_tril).statistic,
            }
        scores = pd.DataFrame(scores).T
        print(scores)
        scores['+gof'] = scores['+fc_corr'] + scores['-fc_diff'] + scores['-fcd_ks']
        # save scores
        scores.to_csv(
            os.path.join(sg.out_dir, "scores.csv")
        )
        # this isn't necessary, but just to be sure, set wIE back to 0
        sg.param_lists["wIE"][:] = 0.0

if __name__ == "__main__":
    run_sim_seed_all()
