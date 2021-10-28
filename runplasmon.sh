#!/bin/bash
#SBATCH -a 1-200
#SBATCH -n 1
#SBATCH -o Plasmon.out-%a
module load julia/1.6.1
export frac=$(echo "print($SLURM_ARRAY_TASK_ID/$SLURM_ARRAY_TASK_COUNT)" | julia)
#echo fraction of highest q value is $frac
julia runplasmon.jl $frac 

