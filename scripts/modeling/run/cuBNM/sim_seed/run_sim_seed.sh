#!/bin/bash
# Expects INPUT_DIR, OUTPUT_DIR and PROJECT_DIR in environment variables
cd "$(dirname "$0")"

dataset=$1
sub=$2
SeedMW=$3

if [[ "$dataset" == "pnc" ]]; then
    cmaes_log_path="${OUTPUT_DIR}/sim/${sub}/ctx_parc-schaefer-100_mean001_thresh-1/6maps_schaefer-100_zscore/cmaes_multimaps_gpu/ctx_parc-schaefer-100_hemi-LR_highpass-013_lowpass-none_exc-inter_G_0.5-4_wee_0.05-0.75_wei_0.05-0.75_wie_0_het-wee-wei_SeedMW-${SeedMW}_SeedSim-410_n-81x210.txt"
    sc_path="${INPUT_DIR}/SC/${sub}/ctx_parc-schaefer-100_mean001_thresh-1_desc-strength.txt"
    fc_tril_path="${INPUT_DIR}/FC/${sub}/ctx_parc-schaefer-100_hemi-LR_highpass-013_lowpass-none_exc-inter_desc-FCtril.txt"
if [[ "$dataset" == "imagen" ]]; then
    exit 1 # not implemented
fi

source $PROJECT_DIR/venv/bin/activate
${PROJECT_DIR}/venv/bin/python \
    ./run_sim_seed.py $dataset $cmaes_log_path $sc_path $fc_tril_path