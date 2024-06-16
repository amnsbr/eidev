#!/bin/bash

n_x=$1
export SES=$2
export SC_LEVEL=${3:-"sub"} # sub will use FU2 SC, other option: "ses"
export PARC=${4:-"schaefer-100"}
export EXC_INTER=${5:-"1"} # from FC
export SC_CONFIG=${6:-"mean001"}
export n_runs=${7:-"2"}
export dataset="imagen"

RUN_DIR=$(dirname "$0")

# specify PROJECT_DIR, INPUT_DIR and OUTPUT_DIR
# note that the following scripts which defines these environment
# variables is not included in the repository
source ./set_env.sh

module load Stages/2024 Python NVHPC GSL/2.7 # specific to JURECA-DC
source ${PROJECT_DIR}/venv/bin/activate

n_subs=$(( $n_x * 12 ))
all_subs=$(python ${RUN_DIR}/main.py get_subs $dataset \
    -n_subs=$n_subs -fc_ses=$SES -sc_level=$SC_LEVEL \
    -n_runs=$n_runs -parc=$PARC -sc_config=$SC_CONFIG \
    -exc_inter=$EXC_INTER)
# convert it to an array
IFS=" " read -r -a all_subs <<< "$all_subs"
# run the jobs (one job (node) per 12 subjects)
for ((i=0;i<$n_x;i++)); do
    start_idx=$(( i * 12 ))
    curr_subs="${all_subs[@]:start_idx:12}"
    curr_job_id=$(sbatch --parsable ${RUN_DIR}/run_cmaes.sbatch $curr_subs)
    echo "$curr_job_id is running $curr_subs"
done
