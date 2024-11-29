#!/bin/bash
# Runs recovery optimizations with the synethetic data as the target
# EXPECTS PROJECT_DIR in the environment variables
# This code was run on the Raven HPC cluster
#SBATCH -J run_recovery
#SBATCH -o ./logs/%j.out
#SBATCH -e ./logs/%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=18
#SBATCH --constraint="gpu"
#SBATCH --gres=gpu:a100:4
#SBATCH --time=24:00:00

RESULTS_DIR="${PROJECT_DIR}/results"
# specify INPUT_DIR and OUTPUT_DIR
# note that the following scripts which defines these environment
# variables is not included in the repository
source ${RESULTS_DIR}/set_env.sh

module load cuda/11.8-nvhpcsdk anaconda/3/2023.03
source ${PROJECT_DIR}/venv/bin/activate

for opt_seed in 1 2; do
    for sim_seed in 410 40; do
        srun --exclusive -n 1 --gres=gpu:1 \
            ${PROJECT_DIR}/venv/bin/python \
            ./run_recovery.py $sim_seed $opt_seed &
    done
done

wait