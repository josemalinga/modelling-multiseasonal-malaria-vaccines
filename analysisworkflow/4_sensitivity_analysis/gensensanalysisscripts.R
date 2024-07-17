########################################
# script gensensanalysisscripts.R
#
# creates scripts for running sensitvity analysis and submits the jobs to cluster 
# INPUTS:
#   exp: experiment name
#   chunk_size: batch size for simulation submission

# OUTPUTS:
#	- OM scenario xml files and simulations in GROUP folder

########################################

gensensanalysisscripts <- function(exp, predicted, scale, manual){
  
  user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
  
  dir.create(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/JOB_OUT"))

  GROUP = "/scicore/home/penny/GROUP/M3TPP/"
  
  SIM_FOLDER=paste0(GROUP,exp,"/")
  ERROR_FOLDER=paste0(SIM_FOLDER,"gp/trained/sensitivity/err/")
  dir.create(ERROR_FOLDER)
  
  file.copy(paste0("/scicore/home/penny/",user,"/M3TPP/analysisworkflow/4_sensitivity_analysis/GP_sens_workflow.sh"), 
            paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_JOBS/GP_sens_workflow.sh"),overwrite=TRUE)

  
  sink(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_JOBS/job_sens_GP.sh"))
  
  cat("#!/bin/bash","\n", sep ="")
  cat("#SBATCH --job-name=sensanalysis","\n", sep ="")
  cat("#SBATCH --account=penny","\n", sep ="")
  #cat("#SBATCH -o ",ERROR_FOLDER,"%A_%a.out","\n", sep ="")
  cat("#SBATCH --output=/dev/null","\n", sep ="")
  cat("#SBATCH --error=/dev/null","\n", sep ="")
  cat("#SBATCH --mem=20G","\n", sep ="")
  cat("#SBATCH --qos=30min","\n", sep ="")
  cat("#SBATCH --cpus-per-task=1","\n", sep ="")
  cat("###########################################","\n", sep ="")
  cat("ml purge","\n", sep ="")
  cat("ml R/3.6.0-foss-2018b","\n", sep ="")
  
  cat("GP_DIR=$1","\n", sep ="")
  cat("PARAM_RANGES_FILE=$2","\n", sep ="")
  cat("SENS_DEST_DIR=$3","\n", sep ="")
  cat("SCALE=$4","\n", sep ="")
  cat("MANUAL=$5","\n", sep ="")
  
  cat("# IMPORTANT: the number of files must equal to the task array length (index starts at 0)","\n", sep ="")
  cat("gp_files=(${GP_DIR}*.RData)","\n", sep ="")
  cat("#i=0","\n", sep ="")
  cat("#echo ${split_files[$i]}","\n", sep ="")
  
  
  cat("# Select scenario file in array","\n", sep ="")
  cat("ID=$(expr ${SLURM_ARRAY_TASK_ID} - 1)","\n", sep ="")
  cat("gp_file=${gp_files[$ID]}","\n", sep ="")
  cat("echo \"Postprocessing for $gp_file\" ","\n", sep ="")
  cat("echo \"Scale arg0: $SCALE\" ","\n", sep ="")
  
  cat("Rscript ../../../analysisworkflow/4_sensitivity_analysis/sens_GP.R $gp_file $PARAM_RANGES_FILE $SENS_DEST_DIR $SCALE $MANUAL","\n", sep ="")
  
  sink()
  
  setwd(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_JOBS/"))
  
  GP_folder = paste0(SIM_FOLDER,"gp/trained/" ,predicted ,"/")
  
  if (manual) {
    param_ranges_file = paste0(SIM_FOLDER,"param_ranges_manual.RData")
  } else {
    param_ranges_file = paste0(SIM_FOLDER,"param_ranges.RData")
  }
  
  sens_folder = paste0(SIM_FOLDER,"gp/trained/sensitivity/")
  
  sys_command = paste("bash GP_sens_workflow.sh", GP_folder, param_ranges_file, sens_folder, scale, manual)
  
  # Run  command
  system(sys_command)

}
