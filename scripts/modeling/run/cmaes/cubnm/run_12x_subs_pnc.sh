#!/bin/bash
n_x=$1
export SC_SUB=${2:-"fc_sub"} # other options: "micamics"
export FC_LEVEL=${3:-"sub"}
export PARC=${4:-"schaefer-100"}
export MAPS_NAME=${5:-"6maps"}
export EXC_INTER=${6:-"1"} # from FC
export SC_CONFIG=${7:-"mean001"}
export n_runs=${8:-"2"}
export subsample=${9:-"false"}
export dataset="pnc"

export RUN_DIR=$(dirname "$0")

# specify PROJECT_DIR, INPUT_DIR and OUTPUT_DIR
# note that the following scripts which defines these environment
# variables is not included in the repository
source $RUN_DIR/set_env.sh

module load cuda/11.8-nvhpcsdk anaconda/3/2023.03
source ${PROJECT_DIR}/venv/bin/activate

get_subs_args="-n_subs=$n_subs -sc_sub=$SC_SUB -fc_level=$FC_LEVEL -n_runs=$n_runs -parc=$PARC -sc_config=$SC_CONFIG -exc_inter=$EXC_INTER -maps_name=$MAPS_NAME"
if [[ "$subsample" == "true" ]]; then
    get_subs_args="$get_subs_args -subsample"
fi

n_subs=$(( $n_x * 12 ))
all_subs=$(${PROJECT_DIR}/venv/bin/python ${RUN_DIR}/main.py get_subs $dataset $get_subs_args)
# convert it to an array
IFS=" " read -r -a all_subs <<< "$all_subs"
# run the jobs (one job (node) per 12 subjects)
for ((i=0;i<$n_x;i++)); do
    start_idx=$(( i * 12 ))
    curr_subs="${all_subs[@]:start_idx:12}"
    curr_job_id=$(sbatch --parsable --export=ALL ${RUN_DIR}/run_cmaes.sbatch $curr_subs)
    echo "$curr_job_id is running $curr_subs"
done
