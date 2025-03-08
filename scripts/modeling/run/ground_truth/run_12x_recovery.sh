#!/bin/bash
n_x=$1
n_seeds=$(( $n_x * 12 ))
cd $(dirname "$0")

# list all the seeds
all_seeds=$(cat $(dirname "$0")/sim_seeds.txt)
IFS=" " read -r -a all_seeds <<< "$all_seeds"

seeds_list=""
counter=0
for seed in "${all_seeds[@]}"; do
    if [[ -d ./recovery/sim-${seed}_opt-1 ]] && [[ -d ./recovery/sim-${seed}_opt-2 ]]; then
        continue
    fi
    seeds_list+=" $seed"
    counter=$(( $counter + 1 ))
    if [[ $counter -ge $n_seeds ]]; then
        break
    fi
done
seeds_list=${seeds_list#" "}
IFS=" " read -r -a seeds_list <<< "$seeds_list"


# run the jobs (one job (node) per 12 seeds)
for ((i=0;i<$n_x;i++)); do
    start_idx=$(( i * 12 ))
    curr_seeds="${seeds_list[@]:start_idx:12}"
    curr_job_id=$(sbatch --parsable --export=ALL $(dirname $0)/run_recovery.sbatch $curr_seeds)
    echo "$curr_job_id is running $curr_seeds"
done
