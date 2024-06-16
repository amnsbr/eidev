#!/bin/bash
#SBATCH -J run_cmaes
#SBATCH -o ./logs/%j.out
#SBATCH -e ./logs/%j.err
#SBATCH -A eidev
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=06:00:00
#SBATCH --partition=dc-gpu
#SBATCH --gres=gpu:4
#SBATCH --mail-user=a.saberi@fz-juelich.de 
#SBATCH --mail-type=FAIL 


subs=${@}

RUN_DIR=$(dirname "$0")

# specify PROJECT_DIR, INPUT_DIR and OUTPUT_DIR
# note that the following scripts which defines these environment
# variables is not included in the repository
source ./set_env.sh

module load Stages/2024 Python NVHPC GSL/2.7 # specific to JURECA-DC
source ${PROJECT_DIR}/venv/bin/activate


for sub in $subs; do
    # specificy FC and SC subject and session
    if [[ "$dataset" == "pnc" ]]; then
        args="$dataset -fc_sub=$sub -sc_sub=$SC_SUB"
    else
        args="$dataset -fc_sub=$sub -fc_ses=$SES -sc_level=$SC_LEVEL"
    fi
    for (( SeedMW = 1; SeedMW <= n_runs; SeedMW++ )); do
        # assign 1 gpu to each job running one optimization run
        srun --exclusive -n 1 --gres=gpu:1 \
            ${RUN_DIR}/venv/bin/python \
            /p/project/eidev/experiments/240422_pnc_imagen_cuBNM/main.py run \
                $args -opt_seed=$SeedMW -parc=$PARC \
                -sc_config=$SC_CONFIG -exc_inter=$EXC_INTER &
    done
done

wait