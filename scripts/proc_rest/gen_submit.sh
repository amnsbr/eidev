#!/bin/bash
# Prints out the .submit instructions for HTCondor jobs given the 
# dataset and the number of new subjects
# Expects PROJECT_DIR, PNC_PROJECT_DIR and IMAGEN_PROJECT_DIR
# in environment variables
# Usage: bash gen_submit.sh <dataset> [optional: <n_subjects> <CPUS> <skip_excluded>] | condor_submit

dataset=$1
n_subjs=$2
CPUS=${3:-'1'}
skip_excluded=${4:-true}

RAM='12G'
DISK='40G'

LOGS_DIR="${PROJECT_DIR}/logs/proc_rest"
if [[ "$dataset" == "pnc" ]]; then
    OUTPUT_DIR="${PNC_PROJECT_DIR}/output"
fi
# Note: for imagen OUTPUT_DIR and FC_DIR will be defined in the for loop depending on $ses
FC_DIR="${OUTPUT_DIR}/FC"
if [[ "$dataset" == "imagen" ]]; then
    FS_DIR="${IMAGEN_PROJECT_DIR}/output/freesurfer"
else
    FS_DIR="${OUTPUT_DIR}/freesurfer"
fi

# run it only for subjects with freesurfer output
all_subjects=$(ls ${FS_DIR} | grep sub)

# run on all if n_subjs not provided
if [ ! $n_subjs ]; then
    n_subjs=$(echo ${all_subjects} | wc -w)
fi


# create the logs dir if it doesn't exist
[ ! -d "${LOGS_DIR}" ] && mkdir -p "${LOGS_DIR}"

# print the .submit header
printf "# The environment
universe       = vanilla
getenv         = True
request_cpus   = ${CPUS}
request_memory = ${RAM}
request_disk   = ${DISK}

# Execution
initial_dir    = ${PROJECT_DIR}
executable     = /bin/bash
\n"

# loop over n_subjs subjects for which the post fmriprep
# does not exist already and is not being processed
current_condor_jobs=$(condor_q --nobatch)
counter=0
for sub_orig in $all_subjects; do
    if [[ "$dataset" == "imagen" ]]; then
        # break sub-xyz_ses-AB to its components 
        # and update dirs accordingly
        IFS="_" read -r sub ses <<< $sub_orig
        ses="${ses#ses-}" # remove ses-
        OUTPUT_DIR="${IMAGEN_PROJECT_DIR}/output/${ses}"
        FC_DIR="${OUTPUT_DIR}/FC"
    else
        sub=$sub_orig
    fi
    if [ ! "$(ls -A ${FC_DIR}/${sub} 2>/dev/null)" ]; then 
        # skip excluded subjects (for other reasons, e.g. fMRI motion)
        # if indicated (default: true)
        if $skip_excluded && [ -f "${OUTPUT_DIR}/excluded/${sub}.txt" ]; then
            continue
        fi
        # Check if freesurfer has been done without errors
        if [ ! -f "${FS_DIR}/${sub}/scripts/recon-all.done" ] || [ -f "${FS_DIR}/${sub}/scripts/recon-all.error" ]; then
            continue
        fi
        # skip subjects that are currently being processed
        if [ $(echo $current_condor_jobs | grep "${PROJECT_DIR}/scripts/proc_rest/proc_rest.sh ${dataset} ${sub_orig}" | wc -l) -gt 0  ]; then
            continue
        fi
        printf "arguments = ${PROJECT_DIR}/scripts/proc_rest/proc_rest.sh ${dataset} ${sub_orig} ${CPUS}\n"
        printf "log       = ${LOGS_DIR}/${dataset}_${sub_orig}_\$(Cluster).\$(Process).log\n"
        printf "output    = ${LOGS_DIR}/${dataset}_${sub_orig}_\$(Cluster).\$(Process).out\n"
        printf "error     = ${LOGS_DIR}/${dataset}_${sub_orig}_\$(Cluster).\$(Process).err\n"
        printf "Queue\n\n"
        counter=$(( $counter + 1 ))
        if [[ $counter -ge $n_subjs ]]; then
            break
        fi
    fi
done