#!/bin/bash
#SBATCH --job-name=create_scenario
#SBATCH --account=smith
#SBATCH --mem=100M
#SBATCH -o /scicore/home/smith/burlyd00/smc_lai/JOB_OUT/run_om.out
#SBATCH --qos=30min
#SBATCH --cpus-per-task=1
###########################################
# Script that creates scenario base xml files
# -o /scicore/home/smith/burlyd00/create_scenario.out
# INPUT:
#	PARAM_TABLE_FILE: text file containing the parameter values for each scenario (1 scenario per row)
#	SCAFFOLD_FILE: .xml scaffold file with @parameter@ locations
#	BASE_FOLDER: directory where a file containing parameter replacement specifications per scenario
#	SCENARIO_FOLDER: directory where a scenario xml file for each row in the paramter table is created
#
# created on 03.04.2019
# monica.golumbeanu@unibas.ch
###########################################
ml purge
ml R/3.6.0-foss-2018b


PARAM_TABLE_FILE=$1
SCAFFOLD_FILE=$2
BASE_FOLDER=$3
SCENARIO_FOLDER=$4

# Extract task ID (corresponds to the line in the parameter table to be used for creating the scenario)
ID=$(expr ${SLURM_ARRAY_TASK_ID})

# Select the parameter names and the scenario  line from the parameter file 
column_names=$(sed -n "1p" < $PARAM_TABLE_FILE)
param_line=$(sed -n "${ID}p" < $PARAM_TABLE_FILE)

# Construct the replacement patterns
Rscript create_scenario.R --column_names $column_names --params $param_line --base_folder $BASE_FOLDER --scenario_folder $SCENARIO_FOLDER --scaffold_file $SCAFFOLD_FILE  
echo "Replacement patterns, base and scenario xml created."


