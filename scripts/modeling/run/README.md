This directory includes the code that was used to perform the main and robustness/supplementary runs of CMA-ES and simulations.

## Using `bnm_cuda`

The majority of simulation-optimizations were performed using the legacy [`bnm_cuda`](https://github.com/amnsbr/bnm_cuda) code on JURECA-DC (an HPC cluster) which uses Slurm.

In `bnm_cuda`, each of the `pnc` and `imagen` subfolders includes the following files:

- `run_12x_subs.sh`: Main script which runs CMA-ES for batches of subjects (12 per job)
- `get_subs.sh`: Determines the next N subjects for whom CMA-ES should be run
- `head.sh`: Determines some constants in addition to input and output paths for a given subject (& session)
- `run_multimaps.sbatch`: Slurm job specification script to run a batch of provided subjects.

Note that some parts of the scripts are specific to JURECA-DC.

These scripts expect `INPUT_DIR`, `BNM_BUILD_PATH`, and `OUTPUT_DIR` in the environment variables:
- `INPUT_DIR` must include:
    - `pnc/SC/$sub`, `pnc/FC/$sub`, `imagen/$ses/SC/$sub`, `imagen/$ses/FC/$sub` and `micamics/SC/group-all` folders created using `scripts/proc_rest` and `scripts/proc_dwi` scripts
    - `${N_MAPS}maps_${PARC}_zscore.txt` created using `scripts/utils/datasets.py`
    - lists of eligible subject IDs (`pnc_subs.txt`, `imagen_subs_FU2.txt`, `imagen_subs_BLnFU2.txt`) and age groups (`pnc_age_groups.txt`)
- `BNM_BUILD_PATH` must include run_CMAES_gpu as compiled in `tools/bnm_cuda/example/run_CMAES.sh`
- The CMAES logs and optimal simulations outputs will be written as .txt files to `OUTPUT_DIR`

Below is a list of commands we used to run the main and robustness/supplementary runs:

### PNC
- Main (subject-specific SC and FC(D) with default configs): `bash ./pnc/run_12x_subs.sh 63 sub sub`
- Effect of SC (Micamics SC and subject-specific FC(D)): `bash ./pnc/run_12x_subs.sh 63 micamics sub`
- Age-groups default (Micamics SC and age-group FC(D)): `bash ./pnc/run_12x_subs.sh 3 micamics age_group`
- Age-groups Schaefer-200 (Micamics SC and age-group FC(D)): `bash ./pnc/run_12x_subs.sh 3 micamics age_group schaefer-200`
- Age-groups 2 maps (Micamics SC and age-group FC(D)): `bash ./pnc/run_12x_subs.sh 3 micamics age_group schaefer-100 2`
- Age-groups 4 maps (Micamics SC and age-group FC(D)): `bash ./pnc/run_12x_subs.sh 3 micamics age_group schaefer-100 4`
- Age-groups including inter-hemispheric in FC(D) (Micamics SC and age-group FC(D)): `bash ./pnc/run_12x_subs.sh 3 micamics age_group schaefer-100 6 false`

### IMAGEN
- Main FU2 (SC: FU2, FC(D): FU2): `bash ./imagen/run_12x_subs.sh 13 FU2 sub`
- Main BL (SC: FU2, FC(D): BL): `bash ./imagen/run_12x_subs.sh 13 BL sub`
- Session-specific BL (SC: BL, FC(D): BL): `bash ./imagen/run_12x_subs.sh 13 BL ses`

## Using `cuBNM`
[`cuBNM`](https://github.com/amnsbr/cuBNM) was developed with a similar core code as [`bnm_cuda`](https://github.com/amnsbr/bnm_cuda) but with an improved interface. We performed three of the additional experiments (`cuBNM/sim_seed`, `cuBNM/delay`, `cuBNM/perturbation`) using cuBNM on Juseless (an HTC cluster) which uses HTCondor. We also provide example code and Slurm job specification scripts for doing the main CMA-ES runs described above using [`cuBNM`](https://github.com/amnsbr/cuBNM) in `cuBNM/cmaes`, though this example code was not used to generate the data reported in the paper and is provided only for the interested reader.

These scripts expect `PROJECT_DIR`, `INPUT_DIR` and `OUTPUT_DIR` in environment variables. `PROJECT_DIR` must include a `./venv` virtual environment with `cuBNM` installed via `pip install git+https://github.com/amnsbr/cuBNM@ab5b8d93479ba9765de0aebf7a24fa198fd769e6`. To run identical simulations as `bnm_cuda` it should be installed after `export CUBNM_NOISE_WHOLE=1`, which enforces precalculation of the entire noise array on CPU rather than noise segments.

### `sim_seed`, `delay` and `perturbation`
Repeat the optimal simulation of each subject with i) 50 alternative random seeds used for Gaussian noise generation, or ii) delayed conduction using 6 alternative velocities, iii) one of the four parameters perturbed by +10% or -10% (this is done for 40 randomly selected subjects). Usage:

1. Generate DAGMan files for optima of PNC subjects using `./cuBNM/<sim_seed|delay|perturbation>/gen_dag.py pnc <path-to-optima-csv>`. This CSV is generated in Figure 2 Jupyter notebook after the main simulations are done (but is not shared in the repository).
2. Run the DAGMan files using `condor_submit_dag <path-to-dag-file>`

### `cmaes`
Runs the main PNC and IMAGEN CMA-ES using `cuBNM` with an equivalent code to the `./bnm_cuda/pnc` and `./bnm_cuda/imagen` main scripts. The simulations (given the same parameters and input data) are identical between `bnm_cuda` and `cuBNM`, but the CMA-ES optimization code used in `cuBNM` is different, and therefore the optimal points resulting from `bnm_cuda` and `cuBNM` will not exactly match.

Example usage:
- PNC main run: `bash ./cuBNM/cmaes/run_12x_subs_pnc.sh 63 fc_sub sub`
- IMAGEN BL main run: `bash ./cuBNM/cmaes/run_12x_subs_imagen.sh 13 BL sub`
- IMAGEN FU2 main run: `bash ./cuBNM/cmaes/run_12x_subs_imagen.sh 13 FU2 sub`