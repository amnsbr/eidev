#!/bin/bash
#SBATCH -J run_cmaes_pnc_cubnm
#SBATCH -o ./logs/%j.out
#SBATCH -e ./logs/%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=18
#SBATCH --constraint="gpu"
#SBATCH --gres=gpu:a100:4
#SBATCH --time=24:00:00 

subs=${@}

export RUN_DIR=$(dirname "$0")

# specify PROJECT_DIR, INPUT_DIR and OUTPUT_DIR
# note that the following scripts which defines these environment
# variables is not included in the repository
source $RUN_DIR/set_env.sh

module load cuda/11.8-nvhpcsdk anaconda/3/2023.03
source ${PROJECT_DIR}/venv/bin/activate

for sub in $subs; do
    # specificy FC and SC subject and session
    if [[ "$dataset" == "pnc" ]]; then
        args="$dataset -fc_subs=$sub -sc_sub=$SC_SUB"
    else
        args="$dataset -fc_subs=$sub -fc_ses=$SES -sc_level=$SC_LEVEL"
    fi
    for (( SeedMW = 1; SeedMW <= n_runs; SeedMW++ )); do
        # assign 1 gpu to each job running one optimization run
        srun --exclusive -n 1 --gres=gpu:1 \
            ${PROJECT_DIR}/venv/bin/python \
            ${RUN_DIR}/main.py run \
                $args -opt_seed=$SeedMW -parc=$PARC \
                -sc_config=$SC_CONFIG -exc_inter=$EXC_INTER \
                -maps_name=$MAPS_NAME &
    done
done

wait