#!/bin/bash

fc_sub=$1
fc_ses=$2

sc_sub=$fc_sub
if [ "$SC_LEVEL" == "ses" ]; then
    sc_ses=$fc_ses
else
    sc_ses="FU2"
fi

EXE="${BNM_BUILD_PATH}/run_CMAES_gpu"

export BNM_CMAES_EARLY_STOP_TOL="0.005"

if [ -z $PARC ]; then
    export PARC="schaefer-100"
fi
if [ -z $SC_CONFIG ]; then
    export SC_CONFIG="mean001"
fi
if [ -z $n_runs_per_subject ]; then
    export n_runs_per_subject=2
fi

maps="${INPUT_DIR}/6maps_${PARC}_zscore.txt"
maps_name=$(basename "$maps" ".txt")

sc_config="${SC_CONFIG}_thresh-1"
sc_path="${INPUT_DIR}/${dataset}/${sc_ses}/SC/${sc_sub}/ctx_parc-${PARC}_${sc_config}_desc-strength.txt"
out_path="${OUTPUT_DIR}/${dataset}/${fc_ses}/${sc_sub}/sc-${sc_ses}_ctx_parc-${PARC}_${sc_config}/${maps_name}"

fc_prefix="ctx_parc-${PARC}_hemi-LR_highpass-013_lowpass-none"
if [[ "$EXC_INTER" == "true" ]]; then
    fc_prefix="${fc_prefix}_exc-inter"
fi
fc_path="${INPUT_DIR}/${dataset}/${fc_ses}/FC/${fc_sub}/${fc_prefix}_desc-FCtril.txt"
fcd_path="${INPUT_DIR}/${dataset}/${fc_ses}/FC/${fc_sub}/${fc_prefix}_desc-FCDtril.txt"

if [[ "$PARC" == "schaefer-100" ]]; then
    regions=100
elif [[ "$PARC" == "schaefer-200" ]]; then
    regions=200
fi

G="0.5-4"
wEE="0.05-0.75"
wEI="0.05-0.75"
wIE="0"
het_params="wee-wei"
dataset="imagen"
TR=2200
fcd_step=2
fcd_window=14
lambda=210
iterMax=81
time_steps=450000
SeedSim=410