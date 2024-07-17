########################################
# script gengridoptimizationscripts.R
#
# creates scripts for running optimization procedure analysis and submits the jobs to cluster 
# INPUTS:
#   exp: experiment name
#   chunk_size: batch size for simulation submission

# OUTPUTS:
#	- database of predictions for many combinations of each continuous parameter
# - database of minimum optimized parameter value need to reach each target outcome

########################################

gengridoptimizationscripts <- function(exp, pred, scale, ngrid, target_range_size){
  
  user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
  
  GROUP <- "/scicore/home/penny/GROUP/M3TPP/"
  
  SIM_FOLDER <- paste0(GROUP, exp, "/")
  ERROR_FOLDER <- paste0(GROUP, exp, "/gp/GP_grid_optimization/err/")
  OPT_FOLDER <- paste0(SIM_FOLDER, "gp/GP_grid_optimization/")
  
  if(!dir.exists(OPT_FOLDER)) dir.create(OPT_FOLDER)
  if(!dir.exists(paste0(OPT_FOLDER, pred, "/"))) dir.create(paste0(OPT_FOLDER, pred, "/"))
  if(!dir.exists(ERROR_FOLDER)) dir.create(ERROR_FOLDER)
  
  file.copy(paste0("/scicore/home/penny/", user, "/M3TPP/analysisworkflow/6_GP_grid_optimization/grid_optimize_workflow.sh"), 
            paste0("/scicore/home/penny/", user, "/M3TPP/Experiments/", exp, "/OM_JOBS/grid_optimize_workflow.sh"), overwrite = TRUE)
  
sink(paste0("/scicore/home/penny/", user, "/M3TPP/Experiments/", exp ,"/OM_JOBS/job_grid_optimize_parameter.sh"))
  
cat("#!/bin/bash", "\n", sep = "")
cat("#SBATCH --job-name=optimization", "\n", sep = "")
cat("#SBATCH --account=penny", "\n", sep = "")
#cat("#SBATCH -o ", ERROR_FOLDER, "%A_%a.out", "\n", sep = "")
cat("#SBATCH --output=/dev/null","\n", sep ="")
cat("#SBATCH --error=/dev/null","\n", sep ="")
cat("#SBATCH --mem=1G","\n", sep ="")
cat("#SBATCH --qos=1day","\n", sep ="")
cat("#SBATCH --cpus-per-task=1","\n", sep ="")
cat("###########################################", "\n", sep = "")

cat("ml purge", "\n", sep = "")
cat("ml R/3.6.0-foss-2018b", "\n", sep = "")

cat("GP_DIR=$1", "\n", sep ="")
cat("SIM_FOLDER=$2", "\n", sep ="")
cat("PRED=$3", "\n", sep = "")
cat("SCALE=$4", "\n", sep = "")
cat("NGRID=$5", "\n", sep = "")
cat("TARGET_RANGE_SIZE=$6", "\n", sep = "")

# For debugging:
# cat("echo $GP_DIR", "\n", sep = "")
# cat("echo $SIM_FOLDER", "\n", sep = "")
# cat("echo $PRED", "\n", sep = "")
# cat("echo $SCALE", "\n", sep = "")
# cat("echo $NGRID", "\n", sep = "")
# cat("echo $TARGET_RANGE_SIZE", "\n", sep = "")

cat("gp_files=(${GP_DIR}*.RData)", "\n", sep = "")
cat("# Select scenario file in array", "\n", sep = "")
cat("ID=$(expr ${SLURM_ARRAY_TASK_ID} - 1)", "\n", sep = "")
cat("gp_file=${gp_files[$ID]}", "\n", sep ="")
cat("echo $gp_file", "\n", sep ="")

cat("ml R/3.6.0-foss-2018b", "\n", sep = "")
cat("Rscript ../../../analysisworkflow/6_GP_grid_optimization/grid_optimize_parameter.R $gp_file $SIM_FOLDER $SCALE $NGRID $TARGET_RANGE_SIZE", "\n", sep = "")

sink()

setwd(paste0("/scicore/home/penny/", user, "/M3TPP/Experiments/", exp, "/OM_JOBS/"))

sys_command = paste("bash grid_optimize_workflow.sh", SIM_FOLDER, pred, scale, ngrid, target_range_size)

# Run  command
system(sys_command)

}
