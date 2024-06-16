#!/bin/bash
# Runs micapipe diffusion processing on the given subject
# + its prerequisite structural processing
# Expects PROJECT_DIR, PNC_PROJECT_DIR, IMAGEN_PROJECT_DIR,
# and SINGULARITY_DIR in environment variables
# Usage: bash proc_dwi.sh <dataset> <participant_label>

dataset=$1
participant_label=$2

# for imagen break down sub-xyz_ses-AB to its components
orig_participant_label=$participant_label
if [[ "$dataset" == "imagen" ]]; then
    IFS="_" read -r participant_label ses <<< $participant_label
    ses="${ses#ses-}" # remove ses-
fi

# set up the directories
MICAPIPE_IMG="${SINGULARITY_DIR}/micapipe-0.1.2.simg"
if [[ "$dataset" == "pnc" ]]; then
    BIDS_DIR="${PNC_PROJECT_DIR}/input/datasets_repo/original/pnc/bids"
    OUTPUT_DIR="${PNC_PROJECT_DIR}/output"
    dwi_main="${BIDS_DIR}/${participant_label}/dwi/${participant_label}_acq-35_dwi.nii.gz,${BIDS_DIR}/${participant_label}/dwi/${participant_label}_acq-36_dwi.nii.gz"
    fs_sub_id=$participant_label
else
    BIDS_DIR="${IMAGEN_PROJECT_DIR}/input/bids/${ses}"
    OUTPUT_DIR="${IMAGEN_PROJECT_DIR}/output/${ses}"
    dwi_main="${BIDS_DIR}/${participant_label}/dwi/${participant_label}_dwi.nii.gz"
    fs_sub_id=${participant_label} # this points to the unarchived freesurfer output
fi

MICAPIPE_DIR="${OUTPUT_DIR}/micapipe"
SC_DIR="${OUTPUT_DIR}/SC"
FS_DIR="${OUTPUT_DIR}/freesurfer"
if [ ! -d "${FS_DIR}/${fs_sub_id}" ]; then
    echo "Run Freesurfer first. Exitting..."
    exit 1
fi

# define shortcut to run the micapipe within
# the singularity image
micapipe () {
    singularity run --cleanenv \
        -B ${PROJECT_DIR},${PNC_PROJECT_DIR},${IMAGEN_PROJECT_DIR} \
        $MICAPIPE_IMG ${*}
}

# run micapipe
if [ "$(ls -A ${SC_DIR}/${participant_label} 2>/dev/null)" ]; then
    echo "micapipe is already done for ${participant_label} of ${dataset}"
else
    echo "running micapipe proc_structural for ${participant_label} ${ses} of ${dataset}"
    micapipe \
        -sub $participant_label \
        -bids $BIDS_DIR \
        -out $OUTPUT_DIR \
        -proc_structural

    echo "running micapipe proc_dwi for ${participant_label} of ${dataset}"
    micapipe \
        -sub $participant_label \
        -bids $BIDS_DIR \
        -out $OUTPUT_DIR \
        -proc_dwi \
        -dwi_main $dwi_main

    echo "running micapipe post_structural for ${participant_label} ${ses} of ${dataset}"
    micapipe \
        -sub $participant_label \
        -bids $BIDS_DIR \
        -out $OUTPUT_DIR \
        -post_structural \
        -atlas schaefer-100,schaefer-200

    echo "running micapipe SC for ${participant_label} ${ses} of ${dataset}"
    micapipe \
        -sub $participant_label \
        -bids $BIDS_DIR \
        -out $OUTPUT_DIR \
        -SC \
        -tracts 10M

    echo "converting SC to labeled CSVs"
    eval "$(conda shell.bash hook)" && \
    conda activate ${PROJECT_DIR}/env && \
    ${PROJECT_DIR}/env/bin/python \
        ${PROJECT_DIR}/scripts/proc_dwi/post_micapipe.py $dataset $orig_participant_label
fi