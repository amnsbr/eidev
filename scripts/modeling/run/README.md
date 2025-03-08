This directory includes the code that was used to perform the main and robustness/supplementary runs of CMA-ES and simulations.

# Running simulations and model optimizations

## Requirements

### Hardware
Nvidia GPUs are required to run the simulations and optimizations. The majority of the analyses in this study was run on Nvidia A100 GPUs, but in some smaller analyses Nvidia GeForce GTX 1080 Ti GPUs were used.

### Software
#### `bnm_cuda`
The legacy [`bnm_cuda`](https://github.com/amnsbr/bnm_cuda) code is needed to perform the main CMA-ES runs. It can be compiled as shown in `tools/bnm_cuda/example/run_CMAES.sh`. See its build dependencies at [here](https://github.com/amnsbr/bnm_cuda?tab=readme-ov-file#build-dependencies). The compiled program `run_CMAES_gpu` should be placed in the `BNM_BUILD_PATH` environment variable.

#### `cubnm`
The [`cubnm`](https://github.com/amnsbr/cubnm) was developed with a similar core code as [`bnm_cuda`](https://github.com/amnsbr/bnm_cuda) but with an improved interface and performance. It should be installed from source via `pip install git+https://github.com/amnsbr/cubnm@v0.0.5.post0`. The requirements for installation from source are listed [here](https://cubnm.readthedocs.io/en/latest/install.html#from-source). Note that to use identical noise arrays as `bnm_cuda` it should be installed after `export CUBNM_NOISE_WHOLE=1`.

## BNM simulation-optimization runs

In `./cmaes` subfolders, we provide the scripts used to run the main and robustness/sensitivity CMA-ES optimization runs for the PNC and IMAGEN datasets.

### Using `bnm_cuda`
The majority of simulation-optimizations were performed using the legacy [`bnm_cuda`](https://github.com/amnsbr/bnm_cuda) code on JURECA-DC and Raven (HPC clusters) which use Slurm.

In `bnm_cuda`, each of the `pnc` and `imagen` subfolders includes the following files:

- `run_12x_subs.sh`: Main script which runs CMA-ES for batches of subjects (12 per job). Usage:
    - `bash ./cmaes/bnm_cuda/pnc/run_12x_subs.sh <n_batches> <sc_level> <fc_level> <parc> <maps_name> <exc_inter> <sc_config> <n_runs> <subsample>`
    - `bash ./cmaes/bnm_cuda/imagen/run_12x_subs.sh <n_batches> <session> <sc_level> <parc> <exc_inter> <sc_config> <n_runs>`
- `get_subs.sh`: Determines the next N subjects for whom CMA-ES should be run
- `head.sh`: Determines some constants in addition to input and output paths for a given subject (& session)
- `run_cmaes.sbatch`: Slurm job specification script to run a batch of provided subjects.

Note that some parts of the scripts are specific to JURECA-DC or Raven HPC systems.

These scripts expect `INPUT_DIR`, `BNM_BUILD_PATH`, and `OUTPUT_DIR` in the environment variables:
- `INPUT_DIR` must include:
    - `pnc/SC/$sub`, `pnc/FC/$sub`, `imagen/$ses/SC/$sub`, `imagen/$ses/FC/$sub` and `micamics/SC/group-all` folders created using `scripts/proc_rest` and `scripts/proc_dwi` scripts
    - `${N_MAPS}maps_${PARC}_zscore.txt` created using `scripts/utils/datasets.py`
    - lists of eligible subject IDs (`pnc_subs.txt`, `pnc_subsample_200.txt`, `imagen_subs_FU2.txt` and `imagen_subs_BLnFU2.txt`)
- The CMAES logs and optimal simulations outputs will be written as .txt files to `OUTPUT_DIR`

### Using `cubnm`
Runs the main PNC and IMAGEN CMA-ES using `cubnm` with an equivalent code to the `./bnm_cuda/pnc` and `./bnm_cuda/imagen` main scripts. 

Usage:
- PNC: `bash ./cmaes/cubnm/run_12x_subs_pnc.sh <n_batches> <sc_level> <fc_level> <parc> <maps_name> <exc_inter> <sc_config> <n_runs> <subsample> <noise_seed>`
- IMAGEN: `bash ./cmaes/cubnm/run_12x_subs_imagen.sh <n_batches> <session> <sc_level> <parc> <exc_inter> <sc_config> <n_runs>`

Note that this code was only used for the following two sensitivity runs on the PNC subsample:
- Schaefer 200, as `cubnm` offers a better performance in running these simulations
- Node-based heterogeneity, as `cubnm` and not `bnm_cuda` supports this feature


## Rerunning optimal simulations with delayed conduction
Using `cubnm` repeat the optimal simulation of each subject with delayed conduction using 6 alternative velocities. Usage: `bash ./delay/run_delay.sh <sub> <cmaes_seed>`. `./delay/gen_dag.py` and `./delay/run_delay.submit` scripts are included to run the simulations on a cluster using HTCondor DAGMan. This script must be run after the main CMA-ES runs are completed.

## Rerunning optimal simulations with different simulation seeds or perturbed parameters
Using `cubnm` repeat the optimal simulation of each subject with i) 50 alternative random seeds used for Gaussian noise generation, or ii) one of the four parameters perturbed by +10% or -10% (this is done for 40 randomly selected subjects). Usage: `python ./sim_seed/run_all.py` and `python ./perturbation/run_all.py`.

## Running the ground truth recovery analysis
First a synthetic ground truth simulation was run using `sbatch ./ground_truth/gen_synthetic.sbatch`. Then the recovery CMA-ES optimizations (using same/different noise seeds listed in `./ground_truth/sim_seeds.txt`) was run using: `bash ./ground_truth/run_12x_recovery.sh`