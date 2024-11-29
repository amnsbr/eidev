"""
Run the simulation for an example PNC subject while saving
the timeseries of model states
"""
import os
import sys
import numpy as np
from cubnm import sim
import set_env

sub = sys.argv[1]
SeedMW = sys.argv[2]
# extract parameters of optimal simulation from `bnm_cuda` output
cmaes_log_path = os.path.join(
    os.environ['PNC_PROJECT_DIR'], 'output', 'sim', sub, 
    'ctx_parc-schaefer-100_mean001_thresh-1',
    '6maps_schaefer-100_zscore','cmaes_multimaps_gpu',
    f'ctx_parc-schaefer-100_hemi-LR_highpass-013_lowpass-none_exc-inter_G_0.5-4_wee_0.05-0.75_wei_0.05-0.75_wie_0_het-wee-wei_SeedMW-{SeedMW}_SeedSim-410_n-81x210.txt'
)
sc_path = os.path.join(os.environ['PNC_PROJECT_DIR'], 'output', 'SC', sub, 'ctx_parc-schaefer-100_mean001_thresh-1_desc-strength.txt')
# get parameters from CMAES log and associated simulation outputs
param_lists = {}
# get G from CMAES log best parameters
with open(cmaes_log_path, 'r') as f:
    lines = f.read().split('\n')
for line in lines:
    if line.startswith('Best'):
        best_params_str = line
if not best_params_str:
    print(f"{cmaes_log_path} is incomplete")
    exit(0)
best_params_str = best_params_str.replace("Best feasible parameters: ","")
best_params = {s.split("=")[0]:float(s.split("=")[1]) for s in best_params_str.split(",")}
param_lists['G'] = np.array([best_params['G']])
# get wEE and wEI arrays from simulation output
sims_out_prefix = cmaes_log_path.replace('/cmaes', '/sims').replace('.txt', '')
param_lists['wEE'] = np.loadtxt(sims_out_prefix + '_w_EE.txt')[np.newaxis,:]
param_lists['wEI'] = np.loadtxt(sims_out_prefix + '_w_EI.txt')[np.newaxis,:]

# run simulation using cubnm
sg = sim.rWWSimGroup(
    duration=450,
    TR=3,
    window_size=10,
    window_step=2,
    sc=sc_path,
    out_dir=os.path.join(cmaes_log_path.replace('.txt', '_cubnm'), 'ts'),
    bw_params="heinzle2016-3T",
    exc_interhemispheric=True,
    rand_seed=410,
    max_fic_trials=10,
    states_ts=True,
    sim_verbose=True,
)
sg.N = 1
sg.param_lists['G'] = param_lists['G']
sg.param_lists['wEE'] = param_lists['wEE']
sg.param_lists['wEI'] = param_lists['wEI']
sg.run()
sg.save()