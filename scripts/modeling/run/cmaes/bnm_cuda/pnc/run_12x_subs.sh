#!/bin/bash

n_x=$1
export SC_LEVEL=${2:-"sub"} # other options: "micamics"
export FC_LEVEL=${3:-"sub"}
export PARC=${4:-"schaefer-100"}
export MAPS_NAME=${5:-"6maps"}
export EXC_INTER=${6:-"true"} # from FC
export SC_CONFIG=${7:-"mean001"}
export n_runs_per_subject=${8:-"2"}
export subsample=${9:-"false"}
export SeedSim=${10:-"410"}

if [[ "$EXC_INTER" != true ]]; then
    export BNM_EXC_INTERHEM=0
fi
echo "BNM_EXC_INTERHEM" $BNM_EXC_INTERHEM

RUN_DIR=$(dirname "$0")
cd $RUN_DIR

source set_env.sh # which sets INPUT_DIR and OUTPUT_DIR
# also export RUN_DIR to be used in run_cmaes.sbatch
export RUN_DIR=$RUN_DIR

n_subs=$(( $n_x * 12 ))
all_subs=$(bash ${RUN_DIR}/get_subs.sh $n_subs)
# convert it to an array
IFS=" " read -r -a all_subs <<< "$all_subs"
# run the jobs (one job (node) per 12 subjects)
for ((i=0;i<$n_x;i++)); do
    start_idx=$(( i * 12 ))
    curr_subs="${all_subs[@]:start_idx:12}"
    curr_job_id=$(sbatch --parsable --export=ALL ${RUN_DIR}/run_cmaes.sbatch $curr_subs)
    echo "$curr_job_id is running $curr_subs"
done
