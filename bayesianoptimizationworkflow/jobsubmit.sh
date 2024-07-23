#!/bin/bash -l

#SBATCH --job-name="rtss_all_bo"
#SBATCH --error=/scicore/home/penny/malinga/bayesianoptimizationworkflow/sub%A_%a.err
#SBATCH --account=penny
#SBATCH --qos=gpu1day
#SBATCH --time=11:45:00
#SBATCH --mem=64GB
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --hint=nomultithread
#SBATCH --gres=gpu:1
#SBATCH --partition=a100

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

export LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

eval "$(/scicore/home/penny/malinga/miniconda3/condabin/conda shell.bash hook)"
conda activate bayesoptim

python main.py