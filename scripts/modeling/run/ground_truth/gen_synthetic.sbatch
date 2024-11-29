#!/bin/bash
# Runs synthetic ground truth simulation
# EXPECTS PROJECT_DIR in the environment variables
# This code was run on the Raven HPC cluster
#SBATCH -o ./logs/%j.out
#SBATCH -e ./logs/%j.err
#SBATCH -D ./
#SBATCH -J run_synthetic
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=gpudev
#SBATCH --time=00:15:00
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000

RESULTS_DIR="${PROJECT_DIR}/results"
# specify INPUT_DIR and OUTPUT_DIR
# note that the following scripts which defines these environment
# variables is not included in the repository
source ${RESULTS_DIR}/set_env.sh

module load cuda/11.8-nvhpcsdk anaconda/3/2023.03

source ${PROJECT_DIR}/venv/bin/activate

${PROJECT_DIR}/venv/bin/python ./gen_synthetic.py