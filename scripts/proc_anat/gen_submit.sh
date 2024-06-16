#!/bin/bash
# Prints out the .submit instructions for HTCondor jobs given the 
# dataset and the number of new subjects
# Expects PROJECT_DIR, PNC_PROJECT_DIR and IMAGEN_PROJECT_DIR
# in environment variables
# Usage: bash gen_submit_freesurfer.sh <dataset> [optional: <n_subjects>] | condor_submit

dataset=$1
n_subjs=$2

CPUS='1'
RAM='12G'
DISK='10G'

LOGS_DIR="${PROJECT_DIR}/logs/run_freesurfer"
if [[ "$dataset" == "pnc" ]]; then
    BIDS_DIR="${PNC_PROJECT_DIR}/input/datasets_repo/original/pnc/bids"
    OUTPUT_DIR="${PNC_PROJECT_DIR}/output"
else
    OUTPUT_DIR="${IMAGEN_PROJECT_DIR}/output"
fi
FS_DIR="${OUTPUT_DIR}/freesurfer"
FS_BAK_DIR="${OUTPUT_DIR}/freesurfer_bak"

if [[ "$dataset" == "imagen" ]]; then
    if [ ! -f "${IMAGEN_PROJECT_DIR}/input/imagen_raw/T1w_ses.txt" ]; then
        eval "$(conda shell.bash hook)" && \
        conda activate ${PROJECT_DIR}/env && \
        ${PROJECT_DIR}/env/bin/python \
            ${PROJECT_DIR}/scripts/manage_data/get_imagen_subjects.py T1w
    fi
    all_subjects=$(cat ${IMAGEN_PROJECT_DIR}/input/imagen_raw/T1w_ses.txt)
    IMAGEN_INPUT_DIR="${IMAGEN_PROJECT_DIR}/input/imagen_raw"
else
    all_subjects=$(ls ${BIDS_DIR} | grep sub)
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
# run on all subjects if n_subjs is not provided
if [ ! $n_subjs ]; then
    n_subjs=$(echo ${all_subjects} | wc -w)
fi

# loop over n_subjs subjects for which the SC output
# does not exist already
counter=0
for sub in $all_subjects; do
    if [ ! -d "${FS_DIR}/$sub" ] && [ ! -d "${FS_BAK_DIR}/$sub" ] || [ -f "${FS_DIR}/$sub/scripts/recon-all.error" ]; then
        printf "arguments = ${PROJECT_DIR}/scripts/proc_anat/run_freesurfer.sh ${dataset} ${sub}\n"
        printf "log       = ${LOGS_DIR}/${dataset}_${sub}_\$(Cluster).\$(Process).log\n"
        printf "output    = ${LOGS_DIR}/${dataset}_${sub}_\$(Cluster).\$(Process).out\n"
        printf "error     = ${LOGS_DIR}/${dataset}_${sub}_\$(Cluster).\$(Process).err\n"
        printf "Queue\n\n"
        counter=$(( $counter + 1 ))
        if [[ $counter -ge $n_subjs ]]; then
            break
        fi
    fi
done