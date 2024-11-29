#!/bin/bash
n_subjects=$1

RUN_DIR=$(dirname "$0")

source ${RUN_DIR}/head.sh

function check_sub(){
    # checks if input exists for the subject but CMAES not done
    sub="$1"
    source ${RUN_DIR}/head.sh $sub $SES
    in_exist=false
    if [ -e "$sc_path" ] && [ -e "$fc_path" ] && [ -e "$fcd_path" ]; then
        in_exist=true
    fi
    out_exist=true
    for (( SeedMW = 1; SeedMW <= n_runs_per_subject; SeedMW++ )); do
        CMAES_out_path="${out_path}/cmaes_multimaps_gpu/\
${fc_prefix}_G_${G}_wee_${wEE}_wei_${wEI}_wie_${wIE}_het-${het_params}\
_SeedMW-${SeedMW}_SeedSim-${SeedSim}_n-${iterMax}x${lambda}.txt"
        if [ ! -e "$CMAES_out_path" ]; then
            out_exist=false
        fi
    done
    excluded=false
    excluded_path="${INPUT_DIR}/${dataset}/excluded/${sub}.txt"
    if [ -e "$excluded_path" ]; then
        excluded=true
    fi
    [ "$in_exist" = true ] && [ "$out_exist" = false ] && [ "$excluded" = false ]
}

# list all the subjects
# the lists are not included in the repository
# and contain a list of subjects separated by new lines
all_subs=()
if [ "$SC_LEVEL" == "ses" ]; then
    while IFS= read -r line; do
    all_subs+=("$line")
    done < ${INPUT_DIR}/imagen_subs_BLnFU2.txt
else
    while IFS= read -r line; do
    all_subs+=("$line")
    done < ${INPUT_DIR}/imagen_subs_FU2.txt
fi

sub_list=""
counter=0
for sub in "${all_subs[@]}"; do
    if check_sub $sub; then
        sub_list+=" $sub"
        counter=$(( $counter + 1 ))
        if [[ $counter -ge $n_subjects ]]; then
            break
        fi
    fi
done
sub_list=${sub_list#" "}
echo $sub_list