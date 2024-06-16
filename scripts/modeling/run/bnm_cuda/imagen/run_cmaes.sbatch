#!/bin/bash
#SBATCH -J run_cmaes_imagen
#SBATCH -o ./logs/%j.out
#SBATCH -e ./logs/%j.err
#SBATCH -A eidev
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=24:00:00
#SBATCH --partition=dc-gpu
#SBATCH --gres=gpu:4
#SBATCH --mail-user=a.saberi@fz-juelich.de 
#SBATCH --mail-type=FAIL

ses=$1
shift
subs=${@}

RUN_DIR=$(dirname "$0")

module load Stages/2023 NVHPC GSL/2.7 # specific to JURECA-DC

for sub in $subs; do
    source ${RUN_DIR}/head.sh $sub $ses
    mkdir -p $out_path
    for (( SeedMW = 1; SeedMW <= n_runs_per_subject; SeedMW++ )); do
        CMAES_out_path="${out_path}/cmaes_multimaps_gpu/\
${fc_prefix}_G_${G}_wee_${wEE}_wei_${wEI}_wie_${wIE}_het-${het_params}\
_SeedMW-${SeedMW}_SeedSim-${SeedSim}_n-${iterMax}x${lambda}.txt"
        if [ -e "$CMAES_out_path" ]; then
            continue
        fi
        srun --exclusive -n 1 --gres=gpu:1 \
            $EXE \
            $sc_path $out_path $fc_path $fcd_path $maps \
            $regions \
            $G $wEE $wEI $wIE $het_params \
            $time_steps $TR $fcd_step $fcd_window $SeedSim \
            $lambda $iterMax $SeedMW &
    done
done

wait