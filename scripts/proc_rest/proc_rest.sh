#!/bin/bash
# Runs fmirprep singularity
# Expects PROJECT_DIR, PNC_PROJECT_DIR, IMAGEN_PROJECT_DIR,
# and SINGULARITY_DIR in environment variables
# Usage: bash proc_dwi.sh <dataset> <participant_label> [<CPUS>]

dataset=$1
participant_label=$2
CPUS=${3:-'1'}

# for imagen break down sub-xyz_ses-AB to its components
if [[ "$dataset" == "imagen" ]]; then
    IFS="_" read -r participant_label ses <<< $participant_label
    ses="${ses#ses-}" # remove ses-
fi

# set up the directories
if [[ "$dataset" == "pnc" ]]; then
    BIDS_DIR="${PNC_PROJECT_DIR}/input/datasets_repo/original/pnc/bids"
    OUTPUT_DIR="${PNC_PROJECT_DIR}/output" 
else
    BIDS_DIR="${IMAGEN_PROJECT_DIR}/input/bids/${ses}"
    OUTPUT_DIR="${IMAGEN_PROJECT_DIR}/output/${ses}"
fi

FS_DIR="${OUTPUT_DIR}/freesurfer"
FMRIPREP_DIR="${OUTPUT_DIR}/fmriprep"
POSTFMRIPREP_DIR="${OUTPUT_DIR}/postfmriprep"
FC_DIR="${OUTPUT_DIR}/FC"
[ ! -d "${FMRIPREP_DIR}" ] && mkdir -p ${FMRIPREP_DIR}
[ ! -d "${POSTFMRIPREP_DIR}" ] && mkdir -p ${POSTFMRIPREP_DIR}
[ ! -d "${FC_DIR}" ] && mkdir -p ${FC_DIR}

if [ ! -d "${FS_DIR}/${participant_label}" ]; then
    echo "Run freesurfer first"
    exit 1
fi


WORK_DIR='/tmp'
FMRIPREP_IMG="${SINGULARITY_DIR}/fmriprep-22.0.0.simg"



# run fmriprep
if [ -f "${FMRIPREP_DIR}/${participant_label}/func/${participant_label}_task-rest_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz" ]; then
    echo "fmriprep is already done for ${participant_label} of ${dataset}"
else
    echo "running fmriprep for ${participant_label} ${ses} of ${dataset}"
        singularity run --cleanenv \
            -B ${WORK_DIR},${PROJECT_DIR},${PNC_PROJECT_DIR},${IMAGEN_PROJECT_DIR},${PROJECT_DIR}/tools/freesurfer_license.txt:/opt/freesurfer/license.txt \
            $FMRIPREP_IMG \
            --random-seed 12345 \
            --skull-strip-fixed-seed \
            --output-spaces MNI152NLin6Asym fsaverage:den-10k \
            --task-id 'rest' \
            --fs-subjects-dir $FS_DIR \
            --skip-bids-validation \
            -w $WORK_DIR --clean-workdir \
            --n_cpus $CPUS \
            $BIDS_DIR $FMRIPREP_DIR participant \
            --participant-label $participant_label
    fi
fi

if [ -d "${FC_DIR}/${participant_label}" ]; then
    echo "postfmriprep already done for ${participant_label} of ${dataset}"
else
    echo "running postfmriprep for ${participant_label} of ${dataset}"
    eval "$(conda shell.bash hook)" && \
    conda activate ${PROJECT_DIR}/env && \
    ${PROJECT_DIR}/env/bin/python \
        ${PROJECT_DIR}/scripts/proc_rest/post_fmriprep.py $dataset $participant_label $ses
fi