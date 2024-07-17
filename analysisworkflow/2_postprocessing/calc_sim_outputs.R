##############################
# calc_sim_outputs.R - script which post-processes the OpenMalaria 
#                      simulation result and calculates necessary outputs
#
# INPUTS:
#       om_results_folder: path of the folder containing the OpenMalaria simulation output files
#       split_file: path of the file containing the scenario parameter specifications for each scenario to be processed
#       dest_dir: path of the folder where output files will be saved
#
# OUTPUTS:
#   In the folder specified at dest_dir, the script will create two files:
#       - one file with results aggregated across seeds (agg_*.*)
#       - one file with results for all seeds (seed_*.*)
#
# created on 04.01.2019
# monica.golumbeanu@unibas.ch
##############################

library(dplyr)    
library(rapportools)


user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
source(paste0("/scicore/home/penny/",user,"/M3TPP/analysisworkflow/2_postprocessing/postprocessing_resources.R"))

##### Main part of script: #####


# Read in command arguments
args = commandArgs(TRUE)
om_results_folder = args[1]
split_file = args[2]
dest_dir = args[3]
dosetime0 = as.integer(args[4])

# Create output file names
split_name = basename(split_file)
dest_table_agg = paste(dest_dir, "agg_", split_name, sep="")
dest_table_seeds = paste(dest_dir, "seeds_", split_name, sep="")

cat("Command arguments:")
print(paste("OM folder:", om_results_folder))
print(paste("Split file:", split_file))
print(paste("Split name:", split_name))
print(paste("Dest dir:", dest_dir))
print(paste("Dest dir agg:", dest_table_agg))
print(paste("Dest dir seeds:", dest_table_seeds))
print(paste("Dosetime:", dosetime0))

# Postprocess the OpenMalaria simulations
cat("Run OM postprocessing function:")
postprocess_OM(results_folder = om_results_folder, param_table_file = split_file,
               final_table_dest = dest_table_agg, final_seed_table_dest = dest_table_seeds, dosetime0)
