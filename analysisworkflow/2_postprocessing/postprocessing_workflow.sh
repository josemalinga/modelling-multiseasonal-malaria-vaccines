#!/bin/bash
#SBATCH --job-name=2b_OM_postprocessing
#SBATCH --account=penny
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --mem=10G
#SBATCH --qos=6hours
#SBATCH --cpus-per-task=1
##############################
# Main script for postprocessing OpenMalaria simulations on the cluster. 
# INPUT:
#       SIM_FOLDER = folder with simulation results
#       PARAM_TAB = file with the simulation parameters
#       OUTPUT_FOLDER = folder with postprocessing results
#       DOSE_TIME = integer representing the survey index to consider for 
#                   evaluating intervention impact - UPDATED DECEMBER 2021
# OUTPUT:
#       The script creates the following folders in the specified OUTPUT_FOLDER:
#               split/ = folder containing scenario parameters for each setting
#               processed/ = folder containing processing results
#
# SYNTHAX: 
#       bash postprocessing_workflow.sh SIM_FOLDER DOSE_TIME
# 
#
# created 14.09.2019
# monica.golumbeanu@unibas.ch
#############################
ml purge
ml R/3.6.0-foss-2018b

SIM_FOLDER=$1
DOSE_TIME=$2

PARAM_TAB=$SIM_FOLDER"param_tab.txt"
PARAM_CAT=$SIM_FOLDER"param_ranges_cat.RData"
OM_FOLDER=$SIM_FOLDER"om/"
OUTPUT_FOLDER=$SIM_FOLDER"postprocessing/"
ERROR_FOLDER=$OUTPUT_FOLDER"err/"
SPLIT_FOLDER=$OUTPUT_FOLDER"split/"

# create the necessary folders
mkdir -p $OUTPUT_FOLDER
mkdir -p $SPLIT_FOLDER
mkdir -p $ERROR_FOLDER

# split the parameter table by setting
Rscript ../../../analysisworkflow/2_postprocessing/split_param_table.R $PARAM_TAB $SPLIT_FOLDER $PARAM_CAT

echo  $SPLIT_FOLDER
# Submit postprocessing array job
split_files=(${SPLIT_FOLDER}*.txt)
NUM=${#split_files[@]}
sbatch -W --array=1-$NUM job_postprocessing.sh $SPLIT_FOLDER $OM_FOLDER $OUTPUT_FOLDER $DOSE_TIME $ERROR_FOLDER
