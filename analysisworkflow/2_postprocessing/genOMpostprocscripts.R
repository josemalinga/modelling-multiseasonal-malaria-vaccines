########################################
# script genGPtrainscripts.R
#
# creates scripts for running OM and submits the jobs to cluster 
# INPUTS:
#   exp: experiment name
#   chunk_size: batch size for simulation submission
#
# OUTPUTS:
#	- OM scenario xml files and simulations in GROUP folder
#
# Updates
# narimane.nekkab@swisstph.ch
# 01.10.2021
#
########################################

genOMpostprocscripts <- function(exp, setting,dosetime0, trialweeks){
  
  user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
  
  if(!dir.exists(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/JOB_OUT"))){
    dir.create(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/JOB_OUT"))
  }
  
  GROUP = "/scicore/home/penny/GROUP/M3TPP/"
  ERROR_FOLDER = paste0(GROUP,exp,"/err/")
  
  SIM_FOLDER=paste0(GROUP,exp,"/")
  
  if (setting=="implementation"){
    
    file.copy(paste0("/scicore/home/penny/",user,"/M3TPP/analysisworkflow/2_postprocessing/postprocessing_workflow.sh"), 
                   paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_JOBS/postprocessing_workflow.sh"),overwrite=TRUE)
    
    sink(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_JOBS/job_postprocessing.sh"))
    
    cat("#!/bin/bash","\n", sep ="")
    cat("#SBATCH --job-name=OM_pp","\n", sep ="")
    cat("#SBATCH --account=penny","\n", sep ="")
    #cat("#SBATCH -e ",ERROR_FOLDER,"%A_%a.err","\n", sep ="") #=/dev/null","\n", sep ="")#
    cat("#SBATCH --error=/dev/null","\n", sep ="")
    cat("#SBATCH --mem=5G","\n", sep ="")
    cat("#SBATCH --qos=6hours","\n", sep ="")
    cat("#SBATCH --cpus-per-task=1","\n", sep ="")
    cat("###########################################","\n", sep ="")
    cat("ml purge","\n", sep ="")
    cat("ml R/3.6.0-foss-2018b","\n", sep ="")
    
    cat("INPUT_DIR=$1","\n", sep ="")
    cat("OM_RESULTS_DIR=$2","\n", sep ="")
    cat("DEST_DIR=$3","\n", sep ="")
    cat("DOSE_TIME=$4","\n", sep ="")

    cat("echo \"Input: $INPUT_DIR\"","\n", sep ="")
    cat("echo \"OM dir: $OM_RESULTS_DIR\"","\n", sep ="")
    cat("echo \"Dest dir: $DEST_DIR\"","\n", sep ="")
    cat("echo \"Dosetime: $DOSE_TIME\"","\n", sep ="")

    # IMPORTANT: the number of files must equal to the task array length (index starts at 0)","\n", sep ="")
    cat("split_files=(${INPUT_DIR}*.txt)","\n", sep ="")
    
    cat("# Select scenario file in array","\n", sep ="")
    cat("ID=$(expr ${SLURM_ARRAY_TASK_ID} - 1)","\n", sep ="")
    cat("split_file=${split_files[$ID]}","\n", sep ="")
    cat("echo \"Postprocessing for $split_file\"","\n", sep ="")
    
    cat("Rscript ../../../analysisworkflow/2_postprocessing/calc_sim_outputs.R $OM_RESULTS_DIR $split_file $DEST_DIR $DOSE_TIME","\n", sep ="")
    
    sink()
    
    setwd(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_JOBS/"))
    
    sys_command = paste("sbatch postprocessing_workflow.sh", SIM_FOLDER ,dosetime0)
    
    # Run  command
    system(sys_command)
    
    }else{
      
      file.copy(paste0("/scicore/home/penny/",user,"/M3TPP/analysisworkflow/2_postprocessing/postprocessing_workflow_trial.sh"), 
                paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_JOBS/postprocessing_workflow_trial.sh"),overwrite=TRUE)
      
      sink(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_JOBS/job_postprocessing_trial.sh"))
      
      cat("#!/bin/bash","\n", sep ="")
      cat("#SBATCH --job-name=OM_pp","\n", sep ="")
      cat("#SBATCH --account=penny","\n", sep ="")
      #cat("#SBATCH -e ",ERROR_FOLDER,"%A_%a.err","\n", sep ="") #=/dev/null","\n", sep ="")#
      cat("#SBATCH --error=/dev/null","\n", sep ="")
      cat("#SBATCH --output=/dev/null","\n", sep ="")
      cat("#SBATCH --mem=5G","\n", sep ="")
      cat("#SBATCH --qos=6hours","\n", sep ="")
      cat("#SBATCH --cpus-per-task=1","\n", sep ="")
      cat("###########################################","\n", sep ="")
      cat("ml purge","\n", sep ="")
      cat("ml R/3.6.0-foss-2018b","\n", sep ="")
      
      
      cat("INPUT_DIR=$1","\n", sep ="")
      cat("OM_RESULTS_DIR=$2","\n", sep ="")
      cat("DEST_DIR=$3","\n", sep ="")
      cat("DOSE_TIME=$4","\n", sep ="")

      # IMPORTANT: the number of files must equal to the task array length (index starts at 0)","\n", sep ="")
      cat("split_files=(${INPUT_DIR}*.txt)","\n", sep ="")
      
      cat("# Select scenario file in array","\n", sep ="")
      cat("ID=$(expr ${SLURM_ARRAY_TASK_ID} - 1)","\n", sep ="")
      cat("split_file=${split_files[$ID]}","\n", sep ="")
      cat("echo \"Postprocessing for $split_file\"","\n", sep ="")
      
      cat("Rscript ../../../analysisworkflow/2_postprocessing/calc_sim_outputs_trial.R $OM_RESULTS_DIR $split_file $DEST_DIR $DOSE_TIME","\n", sep ="")
      
      sink()
      
      setwd(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_JOBS/"))
      
      sys_command = paste("sbatch postprocessing_workflow_trial.sh", SIM_FOLDER , trialweeks[1], trialweeks[2])
      
      # Run  command
      system(sys_command)
      
    }
  
}
 
