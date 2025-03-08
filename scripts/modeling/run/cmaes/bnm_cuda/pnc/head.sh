#!/bin/bash

sub=$1
dataset="pnc"

export BNM_BUILD_PATH="${RUN_DIR}/bnm_build"
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
if [ -z $MAPS_NAME ]; then
    export MAPS_NAME="6maps"
fi


if [[ "$MAPS_NAME" == "homo" ]]; then
    maps="homo"
    maps_name="homo"
    het_params="none"
else
    maps="${INPUT_DIR}/${MAPS_NAME}_${PARC}_zscore.txt"
    maps_name=$(basename "$maps" ".txt")
    het_params="wee-wei"
fi

if [[ "$SC_LEVEL" == "micamics" ]]; then
    sc_path="${INPUT_DIR}/micamics/SC/group-all/ctx_parc-${PARC}_approach-median_${SC_CONFIG}_desc-strength.txt"
    out_path="${OUTPUT_DIR}/${dataset}/group-micamics/ctx_parc-${PARC}_${SC_CONFIG}/${sub}/${maps_name}"
else
    sc_path="${INPUT_DIR}/${dataset}/SC/${sub}/ctx_parc-${PARC}_${SC_CONFIG}_thresh-1_desc-strength.txt"
    out_path="${OUTPUT_DIR}/${dataset}/${sub}/ctx_parc-${PARC}_${SC_CONFIG}_thresh-1/${maps_name}"
fi

fc_prefix="ctx_parc-${PARC}_hemi-LR_highpass-013_lowpass-none"
if [[ "$EXC_INTER" == "true" ]]; then
    fc_prefix="${fc_prefix}_exc-inter"
fi
fc_path="${INPUT_DIR}/${dataset}/FC/${sub}/${fc_prefix}_desc-FCtril.txt"
fcd_path="${INPUT_DIR}/${dataset}/FC/${sub}/${fc_prefix}_desc-FCDtril.txt"

if [[ "$PARC" == "schaefer-100" ]]; then
    regions=100
elif [[ "$PARC" == "schaefer-200" ]]; then
    regions=200
fi

G="0.5-4"
wEE="0.05-0.75"
wEI="0.05-0.75"
wIE="0"
TR=3000
fcd_step=2
fcd_window=10
lambda=210
iterMax=81
time_steps=450000