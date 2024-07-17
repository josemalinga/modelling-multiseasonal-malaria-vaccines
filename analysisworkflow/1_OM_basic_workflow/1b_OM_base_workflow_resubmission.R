###############################################
###############################################
###                                         ###
### STEP 1b: RESUBMISSION OF OM SIMULATIONS ###
###                                         ###
###############################################
###############################################

### -------------------------------------------------------------------------
###
### M3TPP PROJECT:
### Creates scripts for running OM and submits the jobs to cluster
###
### INPUTS:
###   exp: experiment name
###   chunk_size: batch size for simulation submission
###
### OUTPUTS:
###	- OM scenario xml files and simulations in GROUP folder
###
### adapted 10.11.2021
### narimane.nekkab@swisstph.ch
###
### -------------------------------------------------------------------------

##############
### HEADER ###
##############

# Clear environment
rm(list = ls())

# Set seed for replication
set.seed(floor(runif(1)*100))

# Library
library(tgp)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory
setwd(paste0("/scicore/home/penny/",user,"/M3TPP"))

#####################################
### RESUBMISSION OM SIMS FUNCTION ###
#####################################

genOMsimscripts_resubmission <- function(exp, QOS, chunk_size){
  
  user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
  user_dir = paste0("/scicore/home/penny/",user,"/M3TPP")
  
  GROUP = "/scicore/home/penny/GROUP/M3TPP/"
  SIM_FOLDER=paste0(GROUP,exp,"/")
  # ERROR_FOLDER=paste0(GROUP,exp,"/err_new/")
  ERROR_FOLDER=paste0(GROUP,exp,"/err_resubmit/")
  
  dir.create(ERROR_FOLDER)
  
  #Generate script to create scenarios and run OM simulations
  sink(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_JOBS/resubmit_OM.sh"))
  
  cat("#!/bin/bash","\n", sep ="")
  cat("#SBATCH --job-name=iTPP2_OMsim","\n", sep ="")
  cat("#SBATCH --account=penny","\n", sep ="")
  cat("#SBATCH --output=/dev/null","\n", sep ="")
  cat("#SBATCH --error=/dev/null","\n", sep ="")
  cat("#SBATCH --mem=2G","\n", sep ="")
  cat("#SBATCH --qos=",QOS,"\n", sep ="")
  cat("#SBATCH --cpus-per-task=1","\n", sep ="")
  cat("#SBATCH --array=1-$NUM%400","\n", sep ="")
  cat("#SBATCH --exclude=shi121","\n", sep ="")
  cat("#SBATCH --exclude=shi122","\n", sep ="")
  cat("#SBATCH --exclude=shi123","\n", sep ="")
  cat("#SBATCH --exclude=shi129","\n", sep ="")
  
  cat("###########################################","\n", sep ="")
  
  cat("PARAM_TABLE_FILE=$1","\n", sep ="")
  cat("SCAFFOLD_FILE=$2","\n", sep ="")
  cat("BASE_FOLDER=$3","\n", sep ="")
  cat("INPUT_DIR=$4","\n", sep ="")
  cat("DEST_DIR=$5","\n", sep ="")
  
  cat("# Load R", "\n", sep ="")
  cat("ml R/3.6.0-foss-2018b", "\n", sep = "")
  
  cat("ID=$(expr ${SLURM_ARRAY_TASK_ID} - 1)","\n", sep ="")
  cat("# echo \"Debug array ID: \" $ID", "\n", sep ="")
  
  # Generate scenarios
  cat("echo \"Generate replacement patterns and scenario xml\"","\n", sep ="")
  
  cat("# Select the parameter names and the scenario  line from the parameter file", "\n", sep ="")
  cat("column_names=$(sed -n \"1p\" < $PARAM_TABLE_FILE)", "\n", sep ="")
  cat("param_line=$(sed -n \"$(expr ${ID} + 2)p\" < $PARAM_TABLE_FILE)", "\n", sep ="")
  
  cat("echo \"Parameter names: \" $column_names", "\n", sep ="")
  cat("echo \"Parameter values: \" $param_line", "\n", sep ="")
  
  cat("# Construct the replacement patterns", "\n", sep ="")
  cat("Rscript ../../../analysisworkflow/1_OM_basic_workflow/create_scenario.R --column_names $column_names --params $param_line --base_folder $BASE_FOLDER --scenario_folder $INPUT_DIR --scaffold_file $SCAFFOLD_FILE","\n", sep ="")
  cat("echo \"Replacement patterns, base and scenario xml created.\"", "\n", sep ="")
  
  # Load OpenMalaria
  cat("# Load OpenMalaria module and change to folder with resource files","\n", sep ="")
  cat("module purge","\n", sep ="")
  cat("ml OpenMalaria/43.0-iomkl-2019.01","\n", sep ="")
  cat("cd /scicore/home/penny/GROUP/M3TPP/OM_schema43","\n", sep ="")
  
  cat("# ml OpenMalaria/38.0-goolf-1.7.20-Python-2.7.11","\n", sep ="")
  cat("# cd /scicore/home/penny/GROUP/M3TPP/OM_schema38","\n", sep ="")
  cat("# module load OpenMalaria/32-RC3-goolf-1.7.20-Python-2.7.11","\n", sep ="")
  cat("# cd /scicore/home/smith/golmon00/OM_schema32","\n", sep ="")
  
  cat("# IMPORTANT: the number of files must equal to the task array length (index starts at 0)","\n", sep ="")
  
  # Run OpenMalaria simulations
  cat("scenario_file=$(echo $param_line | awk '{print $1}')", "\n", sep ="")
  cat("seed=$(echo $param_line | awk '{print $NF}')", "\n", sep ="")
  cat("# echo \"Debug scenario file: \" $scenario_file", "\n", sep ="")
  cat("# echo \"Debug seed: \" $seed", "\n", sep ="")
  
  cat("scenario_file=${INPUT_DIR}${scenario_file}\"_\"${seed}\".xml\"", "\n", sep = "")
  cat("# echo \"Debug scenario file: \" $scenario_file", "\n", sep ="")
  
  cat("# echo \"Running simulation for $scenario_file\"","\n", sep ="")
  
  cat("# Run OpenMalaria on scenario file","\n", sep ="")
  cat("scenario_name=$(basename \"$scenario_file\" \".xml\")","\n", sep ="")
  cat("output1=$DEST_DIR$scenario_name\"_out.txt\"","\n", sep ="")
  cat("output2=$DEST_DIR$scenario_name\"_cts.txt\"","\n", sep ="")
  cat("echo \"Outputs will be written to $output1 and $output2\"","\n", sep ="")
  cat("openMalaria --scenario $scenario_file --output $output1 --ctsout $output2","\n", sep ="")
  cat("echo \"OpenMalaria simulation ended.\"","\n", sep ="")
  cat("parentdir=\"$(dirname \"$INPUT_DIR\")\"","\n", sep ="")
  
  sink()
  
  
  #################################################################### FOR LARGE EXPERIMENTS
  
  
  ####################### Create new param tables
  
  chunk <- chunk_size
  
  # param_tab_new <- read.table(paste0("/scicore/home/penny/GROUP/M3TPP/",exp,"/param_tab_new.txt"), sep= "\t", header = TRUE, as.is = TRUE, stringsAsFactors = FALSE)
  param_tab_new <- read.table(paste0("/scicore/home/penny/GROUP/M3TPP/",exp,"/param_tab_resubmit.txt"), sep= "\t", header = TRUE, as.is = TRUE, stringsAsFactors = FALSE)
  
  n <- nrow(param_tab_new)
  r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
  split_param_tab <- split(param_tab_new,r)
  num_tab <- max(r)
  
  for(j in 0:(num_tab-1)){ 
    
    # table_file = paste0("/scicore/home/penny/GROUP/M3TPP/",exp,"/param_tab_new_",j,".txt")
    table_file = paste0("/scicore/home/penny/GROUP/M3TPP/",exp,"/param_tab_resubmit_",j,".txt")
    
    # Write table to specified destination file
    write.table(split_param_tab[[j+1]], table_file, sep = "\t", quote = FALSE, col.names = TRUE,
                row.names = FALSE)
  }
  
  ####################### Create jobs
  
  # Extract the number of lines in the parameter table
  # no.commands=as.numeric(system(paste0("wc -l < /scicore/home/penny/GROUP/M3TPP/",exp, "/param_tab_new.txt"), intern = TRUE))
  no.commands=as.numeric(system(paste0("wc -l < /scicore/home/penny/GROUP/M3TPP/",exp, "/param_tab_resubmit.txt"), intern = TRUE))
  
  no.bats = no.commands %/% chunk_size
  
  if(no.commands %% chunk_size >0){no.bats = no.bats+1}
  for(j in 0:(no.bats-1)){ #"split" counting automatically starts at 0, so myst start counting from 0 here
    
    sink(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_JOBS/resubmission","_",exp,"_",sprintf("%02i",j),".sh"))
    
    cat("#!/bin/bash\n")
    cat("#SBATCH --job-name=iTPP2_OMsim_resub",sprintf("%02i",j),"\n", sep ="")
    cat("#SBATCH --output=/dev/null","\n", sep ="")
    cat("#SBATCH --error=/dev/null","\n", sep ="")
    cat("#SBATCH --qos=",QOS,"\n", sep ="") 
    cat("PARAM_TABLE_FILE=",SIM_FOLDER,"param_tab_resubmit_",j,".txt","\n\n", sep ="")
    cat("SCAFFOLD_FILE=",SIM_FOLDER,"scaffold.xml","\n", sep ="")
    cat("BASE_FOLDER=",SIM_FOLDER,"base_resubmit_",j,"/","\n", sep ="")
    cat("SCENARIOS_FOLDER=",SIM_FOLDER,"scenarios_resubmit_",j,"/","\n", sep ="")
    cat("OM_FOLDER=",SIM_FOLDER,"om/","\n", sep ="")
    
    cat("NUM=$(wc -l < $PARAM_TABLE_FILE)","\n", sep ="")
    
    cat(" echo \"Creating $NUM-1 scenarios and running OM simulations... \" ","\n", sep ="")
    
    cat("sbatch -W --array=1-$NUM ./run_OM.sh $PARAM_TABLE_FILE $SCAFFOLD_FILE $BASE_FOLDER $SCENARIOS_FOLDER $OM_FOLDER","\n\n", sep ="")
    
    # Close the sink!
    sink()
  }
  
  # either submit run to cluster here or from terminal 
  setwd(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_JOBS/"))
  
  for (j in 0:(no.bats-1)){ 
    sys_command = paste0("sbatch", " resubmission_",exp,"_",sprintf("%02i",j),".sh")
    
    # Run  command
    system(sys_command)
  }
  
} 

##################
### EXPERIMENT ###
##################

# Insert experiment name here
exp ="..."

# Optimize file size
chunk_size = 90000

QOS = "30min"

###########################################################
### GENERATE SCENARIOS AND RUN OPEN MALARIA SIMULATIONS ###
###########################################################
getwd()
# Run
genOMsimscripts_resubmission(exp, QOS, chunk_size)
