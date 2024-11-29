#!/bin/bash

n_x=$1
export SES=$2
export SC_LEVEL=${3:-"sub"} # sub will use FU2 SC, other option: "ses"
export PARC=${4:-"schaefer-100"}
export EXC_INTER=${5:-"1"} # from FC
export SC_CONFIG=${6:-"mean001"}
export n_runs=${7:-"2"}
export dataset="imagen"
export MAPS_NAME="6maps"

export RUN_DIR=$(dirname "$0")

# specify PROJECT_DIR, INPUT_DIR and OUTPUT_DIR
# note that the following scripts which defines these environment
# variables is not included in the repository
source $RUN_DIR/set_env.sh

module load cuda/11.8-nvhpcsdk anaconda/3/2023.03
source ${PROJECT_DIR}/venv/bin/activate

n_subs=$(( $n_x * 12 ))
all_subs=$(${PROJECT_DIR}/venv/bin/python ${RUN_DIR}/main.py get_subs $dataset \
    -n_subs=$n_subs -fc_ses=$SES -sc_level=$SC_LEVEL \
    -n_runs=$n_runs -parc=$PARC -sc_config=$SC_CONFIG \
    -exc_inter=$EXC_INTER -maps_name=$MAPS_NAME)
# convert it to an array
IFS=" " read -r -a all_subs <<< "$all_subs"
# run the jobs (one job (node) per 12 subjects)
for ((i=0;i<$n_x;i++)); do
    start_idx=$(( i * 12 ))
    curr_subs="${all_subs[@]:start_idx:12}"
    curr_job_id=$(sbatch --parsable --export=ALL ${RUN_DIR}/run_cmaes.sbatch $curr_subs)
    echo "$curr_job_id is running $curr_subs"
done
