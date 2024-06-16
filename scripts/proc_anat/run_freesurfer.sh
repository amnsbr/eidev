#!/bin/bash
# Runs freesurfer on a given subject
# Expects PROJECT_DIR, PNC_PROJECT_DIR, IMAGEN_PROJECT_DIR,
# and SINGULARITY_DIR in environment variables
# Usage: bash run_freesurfer.sh <dataset> <subject_id>
dataset=$1
fs_subject_id=$2 # for pnc this is participant_label; for imagen is sub-{}_ses-{}

# specifiy the location of input T1w and the output (SUBJECTS_DIR)
if [[ "$dataset" == "imagen" ]]; then
    export SUBJECTS_DIR="${IMAGEN_PROJECT_DIR}/output/freesurfer"
    t1w_in="${IMAGEN_PROJECT_DIR}/input/imagen_raw/${fs_subject_id}/T1w.nii.gz"
    if [ ! -f "${t1w_in}" ]; then
        # download the input data
        # this is temporary and a special case for imagen to enable
        # running its structural pipeline separately without having to
        # download all the imagen raw files on storage
        eval "$(conda shell.bash hook)" && \
        conda activate ${PROJECT_DIR}/env && \
        ${PROJECT_DIR}/env/bin/python \
            ${PROJECT_DIR}/scripts/manage_data/download_imagen.py $fs_subject_id
    fi
    if [ ! -f "${t1w_in}" ]; then
        # if it still does not exist this subject probably doesn't have a T1w image
        # or the download has failed 
        exit 1
    fi
else
    export SUBJECTS_DIR="${PNC_PROJECT_DIR}/output/freesurfer"
    t1w_in="${PNC_PROJECT_DIR}/input/datasets_repo/original/pnc/bids/${fs_subject_id}/anat/${fs_subject_id}_T1w.nii.gz"
    if [ ! -f "${t1w_in}" ]; then
        # download the input data
        bash ${PROJECT_DIR}/scripts/manage_data/download_input.sh pnc $fs_subject_id
    fi
fi

# create subjects dir if it doesn't exist
[ ! -d $SUBJECTS_DIR ] && mkdir -p $SUBJECTS_DIR

# check if FS is already done without errors
if [ -d "$SUBJECTS_DIR/$fs_subject_id" ] && [ ! -f "$SUBJECTS_DIR/$fs_subject_id/scripts/recon-all.error" ]; then
    echo "Freesurfer already done for ${fs_subject_id}"
    exit 0
fi

# wrapper to run recon-all from singularity
recon-all () {
    singularity exec --cleanenv \
        -B ${PROJECT_DIR},${PNC_PROJECT_DIR},${IMAGEN_PROJECT_DIR},${PROJECT_DIR}/tools/freesurfer_license.txt:/usr/local/freesurfer/license.txt \
        "${SINGULARITY_DIR}/freesurfer-7.1.1.simg" \
        /bin/bash -c \
        "export SUBJECTS_DIR=${SUBJECTS_DIR} && \
        recon-all ${*}"
}

recon-all \
    -subjid $fs_subject_id \
    -all \
    -i $t1w_in \
    -sd $SUBJECTS_DIR