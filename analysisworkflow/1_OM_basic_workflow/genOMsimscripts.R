########################################
# script genOMsimscripts.R
#
# creates scripts for running OM and submits the jobs to cluster 
# INPUTS:
#   exp: experiment name
#   chunk_size: batch size for simulation submission

# OUTPUTS:
#	- OM scenario xml files and simulations in GROUP folder

########################################

genOMsimscripts <- function(exp, chunk_size){
  
  user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
  user_dir = paste0("/scicore/home/penny/",user,"/M3TPP")
  
  dir.create(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/JOB_OUT"))

  GROUP = "/scicore/home/penny/GROUP/M3TPP/"
  
  SIM_FOLDER=paste0(GROUP,exp,"/")
  ERROR_FOLDER=paste0(GROUP,exp,"/err/")
  
  dir.create(SIM_FOLDER)
  dir.create(ERROR_FOLDER)
  
  file.copy(paste0(user_dir,"/Experiments/",exp,"/OM_JOBS/scaffold.xml"), paste0(GROUP,exp,"/scaffold.xml"),overwrite=TRUE)
  file.copy(paste0(user_dir,"/Experiments/",exp,"/OM_JOBS/param_tab.txt"), paste0(GROUP,exp,"/param_tab.txt"),overwrite=TRUE)
  
  #Generate script to create scenarios and run OM simulations
  sink(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_JOBS/run_OM.sh"))
  
cat("#!/bin/bash","\n", sep ="")
cat("#SBATCH --job-name=CT","\n", sep ="")
cat("#SBATCH --account=penny","\n", sep ="")
cat("#SBATCH -o ",ERROR_FOLDER,"%A_%a.out","\n", sep ="")
# cat("#SBATCH -e /scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/JOB_OUT/OM_error.err","\n", sep ="")
cat("#SBATCH --mem=2G","\n", sep ="")
cat("#SBATCH --qos=30min","\n", sep ="")
cat("#SBATCH --cpus-per-task=1","\n", sep ="")
cat("###########################################","\n", sep ="")
cat("# %A_%a.out","\n", sep ="")
cat("# Script for submitting OpenMalaria simulations as jobs to the cluster","\n", sep ="")
cat("# Arguments:","\n", sep ="")
cat("#INPUT_DIR: directory containing the scenario files (.xml)","\n", sep ="")
cat("#DEST_DIR: directory where two output files will be created per OpenMalaria simulation","\n", sep ="")
cat("# Calling the script:","\n", sep ="")
cat("#sbatch --array=1-NUM run_OM.sh PARAM_TABLE_FILE SCAFFOLD_FILE BASE_FOLDER INPUT_DIR DEST_DIR","\n", sep ="")
cat("#","\n", sep ="")
cat("# created on 01.02.2019","\n", sep ="")
cat("# monica.golumbeanu@unibas.ch","\n", sep ="")
cat("###########################################","\n", sep ="")

cat("PARAM_TABLE_FILE=$1","\n", sep ="")
cat("SCAFFOLD_FILE=$2","\n", sep ="")
cat("BASE_FOLDER=$3","\n", sep ="")
cat("INPUT_DIR=$4","\n", sep ="")
cat("DEST_DIR=$5","\n", sep ="")

cat("# Load R", "\n", sep ="")
cat("ml R/3.6.0-foss-2018b", "\n", sep = "")
  
cat("# create destination directory","\n", sep ="")
cat("mkdir -p $DEST_DIR","\n", sep ="")

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
# cat("# Remove error files","\n", sep ="")
# cat(paste0("del \"",ERROR_FOLDER,"%A_%a.out\"","\n"))

sink()


  #generate script for generating scenarios
  # Note that this script is now redundant. Retained here for legacy purposes
  
  # sink(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_JOBS/job_create_scenarios.sh"))
  # 
  # cat("#!/bin/bash","\n", sep ="")
  # cat("#SBATCH --job-name=create_scenarios","\n", sep ="")
  # cat("#SBATCH --account=penny","\n", sep ="")
  # cat("#SBATCH -o /scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/JOB_OUT/create_scenario.out","\n", sep ="")
  # cat("#SBATCH --mem=100M","\n", sep ="")
  # cat("#SBATCH --qos=30min","\n", sep ="")
  # cat("#SBATCH --cpus-per-task=1","\n", sep ="")
  # cat("###########################################","\n", sep ="")
  # cat("# Script that creates scenario base xml files","\n", sep ="")
  # cat("# ","\n", sep ="")
  # cat("# INPUT:","\n", sep ="")
  # cat("#	PARAM_TABLE_FILE: text file containing the parameter values for each scenario (1 scenario per row)","\n", sep ="")
  # cat("#	SCAFFOLD_FILE: .xml scaffold file with @parameter@ locations","\n", sep ="")
  # cat("#	BASE_FOLDER: directory where a file containing parameter replacement specifications per scenario","\n", sep ="")
  # cat("#	SCENARIO_FOLDER: directory where a scenario xml file for each row in the paramter table is created","\n", sep ="")
  # cat("# created on 03.04.2019","\n", sep ="")
  # cat("# monica.golumbeanu@unibas.ch","\n", sep ="")
  # cat("###########################################","\n", sep ="")
  # cat("ml purge","\n", sep ="")
  # cat("ml R/3.6.0-foss-2018b","\n", sep ="")
  # 
  # 
  # cat("PARAM_TABLE_FILE=$1","\n", sep ="")
  # cat("SCAFFOLD_FILE=$2","\n", sep ="")
  # cat("BASE_FOLDER=$3","\n", sep ="")
  # cat("SCENARIO_FOLDER=$4","\n", sep ="")
  # 
  # cat("# Extract the number of lines in the parameter table","\n", sep ="")
  # cat("NUM=$(wc -l < $PARAM_TABLE_FILE)","\n", sep ="")
  # 
  # # Create the scenarios
  # cat("echo \"Creating $NUM-1 scenarios ...\" ","\n", sep ="")
  # cat("for (( ID=2; ID<=$NUM; ID++ ))","\n", sep ="")
  # cat("  do  ","\n", sep ="")
  # cat("# Select the parameter names and the scenario  line from the parameter file ","\n", sep ="")
  # cat("column_names=$(sed -n \"1p\" < $PARAM_TABLE_FILE)","\n", sep ="")
  # cat("param_line=$(sed -n \"${ID}p\" < $PARAM_TABLE_FILE)","\n", sep ="")
  # 
  # cat("# Construct the replacement patterns","\n", sep ="")
  # cat("Rscript ../../../analysisworkflow/1_OM_basic_workflow/create_scenario.R --column_names $column_names --params $param_line --base_folder $BASE_FOLDER --scenario_folder $SCENARIO_FOLDER --scaffold_file $SCAFFOLD_FILE","\n", sep ="")  
  # cat("echo \"Replacement patterns, base and scenario xml created.\"","\n", sep ="")
  # cat("done","\n", sep ="")
  # 
  # cat("# Cleanup","\n", sep ="")
  # cat("rm -r $BASE_FOLDER","\n", sep ="")
  # 
  # sink()
  
  
  
  
  
  # Extract the number of lines in the parameter table
  no.commands=as.numeric(system(paste0("wc -l < /scicore/home/penny/GROUP/M3TPP/",exp, "/param_tab.txt"), intern = TRUE))
  
  no.bats = no.commands %/% chunk_size
  if(no.commands %% chunk_size >0){no.bats = no.bats+1}
  for(j in 0:(no.bats-1)){ #"split" counting automatically starts at 0, so myst start counting from 0 here
    sink(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_JOBS/submission","_",exp,"_",sprintf("%02i",j),".sh"))
    cat("#!/bin/bash\n")
    cat("#SBATCH -o /scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/JOB_OUT/submission",sprintf("%02i",j),".out","\n",sep ="")
    
    cat("#SBATCH --qos=",QOS,"\n", sep ="") #uncomment to change the queue from 6h to 30min (suitable for small jo)
    
    cat("PARAM_TABLE_FILE=",SIM_FOLDER,"param_tab_",j,".txt","\n\n", sep ="")
    cat("SCAFFOLD_FILE=",SIM_FOLDER,"scaffold.xml","\n", sep ="")
    cat("BASE_FOLDER=",SIM_FOLDER,"base_",j,"/","\n", sep ="")
    cat("SCENARIOS_FOLDER=",SIM_FOLDER,"scenarios_",j,"/","\n", sep ="")
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
  sys_command = paste0("sbatch", " submission_",exp,"_",sprintf("%02i",j),".sh")
  
  # Run  command
  system(sys_command)
 }
  
} 
