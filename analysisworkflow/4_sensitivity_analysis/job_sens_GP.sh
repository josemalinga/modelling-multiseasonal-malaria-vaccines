#!/bin/bash
#SBATCH --job-name=sens_GP
#SBATCH --account=smith
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --mem=20G
#SBATCH --qos=30min
#SBATCH --cpus-per-task=1

  cat("#SBATCH --output=/dev/null","\n", sep ="")
  cat("#SBATCH --error=/dev/null","\n", sep ="")

#######################################
# script for sensitivity analysis
#
# created 17.06.2019
# monica.golumbeanu@unibas.ch
######################################
ml purge
ml R/3.6.0-foss-2018b

GP_DIR=$1
PARAM_RANGES_FILE=$2
SENS_DEST_DIR=$3

# IMPORTANT: the number of files must equal to the task array length (index starts at 0)
gp_files=(${GP_DIR}*.RData)
#i=0
#echo ${split_files[$i]}


# Select scenario file in array
ID=$(expr ${SLURM_ARRAY_TASK_ID} - 1)
gp_file=${gp_files[$ID]}
echo "Postprocessing for $gp_file"

Rscript sens_GP_pppy_lowHL.R $gp_file $PARAM_RANGES_FILE $SENS_DEST_DIR
