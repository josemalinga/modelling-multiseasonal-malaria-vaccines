#############################################################################
#############################################################################
#######
####### PROBLEM FUNCTION
#######
####### Adapting OpenMalaria model calibration approach
####### to mAb Phase 2 clinical trial in silico modelling
####### to identify EIR, diagnostic sensitivity, and access (for control arm)
####### and initial efficacy, half-life, and decay (for treatment arms)
#######
#######
#############################################################################
#############################################################################

##############
### HEADER ###
##############

# Close connections
closeAllConnections()

# Clear environment
rm(list = ls())

# Set seed for replication
set.seed(42)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory (set to work in SciCore home directory)
wd = paste0("/scicore/home/penny/",user,"/bayesianoptimizationworkflow/")
setwd(wd)
source(paste0(wd,"clinicaltrial_workflow/functions/convert_access.R"))

# Libraries
library(data.table)
library(stringr)
library(dplyr)

##################
### EXPERIMENT ###
##################

# Arms to simulate (names need to correspond to all data and files)
#arms = c("control")
arms = c("RTSS_placebo", "RTSS_trial","control") # in order to fit to all the arms using all the data

# one country at a time ("BF", "Mali")
country = season = "Mali"

# Insert experiment name here
experiment = exp = paste0("rtss_bo_validation_",country)

######################
### CREATE FOLDERS ###
######################

# Create experiment files
FILES = paste0(wd,"Experiments_FILES/",exp)

# Folder directory for all simulation outputs
OUTPUTS = paste0(wd,"Experiments_OUTPUTS/",exp)

# Create folders
if(!dir.exists(OUTPUTS)){dir.create(OUTPUTS)}
if(!dir.exists(FILES)){dir.create(FILES)}
if(!dir.exists(paste0(OUTPUTS,"/JOBS"))){dir.create(paste0(OUTPUTS,"/JOBS"))}

##################
### PARAMETERS ###
##################

# Load X values
x_parameters = read.csv(paste0(wd, "/X.txt"), sep = " ", header=F)
# write.table(x_parameters, "X.txt", row.names = F, col.names = F)

# Dictionary (dummy trial)
colnames(x_parameters) = c("efficacy4","efficacy5", "halfsmc", "kdecay")

# includes boosting efficacy4, efficacy5 and halflife of SMC in parameter list
x_param_range_OM = x_parameters

x_param_range_OM$efficacy = 0.91
x_param_range_OM$halflife = 223

# include EIr in parameter list
if (country=="Mali"){
x_param_range_OM = data.frame(tidyr::crossing(x_param_range_OM, EIR = seq(from=10, to=22, by=4)))
}

if (country=="BF"){
  x_param_range_OM = data.frame(tidyr::crossing(x_param_range_OM, EIR = seq(from=48, to=96, by=12)))
}

# include access level in parameter list
x_param_range_OM = data.frame(tidyr::crossing(x_param_range_OM, access = c(0.45)))

# Convert access to 5 day probabilities
x_param_range_OM$access = round(convert_access(x_param_range_OM$access * 100), digits = 4)
x_param_range_OM$access = ifelse(x_param_range_OM$access < 0.04, 0.04, x_param_range_OM$access)

#Specify parameter ranges
x_param_range_OM$efficacy4 = scales::rescale(x_param_range_OM$efficacy4, to = c(0.5,0.9))
x_param_range_OM$efficacy5 = scales::rescale(x_param_range_OM$efficacy5, to = c(0.4,0.8))
x_param_range_OM$halfsmc = scales::rescale(x_param_range_OM$halfsmc, to = c(20,32))
x_param_range_OM$kdecay = scales::rescale(x_param_range_OM$kdecay, to = c(2,6))

# using the vaccine weibull decay for bounds
L= 7.3; k = 0.69; eff_init = 91.1
eff_weibull= eff_init * exp( -(12/L)^k * log(2) )

x_param_range_OM = x_param_range_OM[x_param_range_OM$efficacy4<=eff_weibull,]
x_param_range_OM = x_param_range_OM[x_param_range_OM$efficacy5<=x_param_range_OM$efficacy4,]

x_param_range_OM

x_parameters = bind_rows(replicate(10,x_parameters, simplify = FALSE))
x_parameters = x_parameters[c(1:nrow(x_param_range_OM)),]

x_parameters = data.frame(x_parameters)
                          
########################
### PROBLEM FUNCTION ###
########################

# Create function
myFunction = function(x_parameters){
  
  #############################
  ### PARAMETER TRANSLATION ###
  
  # Convert
  x_parameters_OM = x_param_range_OM[,c("efficacy4","efficacy5","halfsmc", "efficacy",
                                        "halflife","EIR","access", "kdecay")]
  x_parameters_OM 
   
  #################################
  ### GENERATE PARAMETER TABLES ###
  
  # Number of seeds
  nseeds = 10
  
  for(i in arms){
  # Source functions
  source(paste0(wd,"clinicaltrial_workflow/functions/getOMparams_",i,".R"))
  source(paste0(wd,"clinicaltrial_workflow/functions/convert_access.R"))

  # Generate parameter table for OM sims across all arms
    getOMparams(wd, i, x_parameters_OM, nseeds, season)
   }
  
  ########################################
  ### SUBMIT OPEN MALARIA SIMULATIONS ###
  
  # Submit OM simulations
  for(i in arms){
    
    # Generate sbatch command with wait time
    sys_command = paste0("#!/bin/bash
      set -e
      date
      
      sbatch -W --job-name=rtss_bo_OM_sims --qos=30min --output=",
                         OUTPUTS,"/JOBS/",i,"_%a.out --error=",
                         OUTPUTS,"/JOBS/",i,"_%a.err --array=1-",
                         nrow(x_parameters)*nseeds," ",
                         wd,"clinicaltrial_workflow/functions/run_OM.sh ",
                         OUTPUTS,"/param_table_",i,".txt ",
                         wd,"clinicaltrial_workflow/data/scaffold_",i,".xml ",
                         OUTPUTS,"/base_",i,"/ ",
                         OUTPUTS,"/scenarios_",i,"/ ",
                         OUTPUTS,"/om_",i,"/;",
                         "
      wait
      
      date
      echo \"all jobs finished\"")
    
    # Run command
    system(sys_command, wait = T)
    
  }
  
  # Remove scenario files when completed
  for(i in arms){
    scenario_files = list.files(paste0(OUTPUTS,"/scenarios_",i,"/"), full.names = T)
    #unlink(scenario_files)
    base_files = list.files(paste0(OUTPUTS,"/base_",i,"/"), full.names = T)
    #unlink(base_files)
  }
  
  ############################
  ### RUN RSS CALCULATIONS ###
  
  # Submit RSS jobs
  for(i in arms){
    
    # Generate sbatch command with wait time
    sys_command = paste0("#!/bin/bash
      set -e
      date
      
      sbatch -W --job-name=rtss_bo_RSS --qos=30min --output=",
                         OUTPUTS,"/JOBS/",i,"_%a.out --error=",
                         OUTPUTS,"/JOBS/",i,"_%a.err --array=1-",
                         nrow(x_parameters)*nseeds," ",
                         wd,"clinicaltrial_workflow/functions/run_RSS.sh ",
                         OUTPUTS,"/param_table_",i,".txt ",
                         OUTPUTS,"/om_",i,"/ ",
                         OUTPUTS,"/RSS_",i,"/ ",
                         i," ",
                         wd,paste0("clinicaltrial_workflow/data/DUMMY_RTSS_",state,"_12MONTHS.RData "),
                         "6;",# follow-up duration in months
                         "
      wait
      
      date
      echo \"all jobs finished\"")
    
    # Run command
    system(sys_command, wait = T)
    
  }
  
  # When finished, load results
  RSS = array(dim = c(nrow(x_parameters),length(arms)))
  # RSS = list()
  for(i in arms){
    # Get file names and load
    file_list = list.files(paste0(OUTPUTS,"/RSS_",i), full.names = T)
    data = lapply(file_list, data.table::fread)
    data = rbindlist(data)
    # Get scenario number
    data$id = str_replace(data$Scenario_Name,paste0(exp,"_"),"")
    id_list = str_split(data$id, "_", n=2)
    data$order=as.numeric(unlist(lapply(id_list, `[[`, 1)))
    # Average by seeds and order
    data = data %>%
      dplyr::group_by(order) %>% 
      dplyr::summarise(RSS = mean(RSS)) %>% 
      dplyr::ungroup() %>% 
      dplyr::arrange(order)
    
    # If missing data, add NA
    data_all = data.frame(order = seq(1:nrow(x_parameters)))
    data_all = merge(data_all, data, by = "order", all.x = T)
    
    # replace all missing data with -1000 to allow the code to run
    data_all[is.na(data_all)] = -1000
    
    # Save RSS
    RSS[,which(i == arms)] = data_all$RSS
    
    # Remove files
    #unlink(file_list)
  }
  
  # Save results
  results = cbind(x_parameters_OM,RSS)
  colnames(results)[(ncol(x_parameters_OM)+1):ncol(results)] = paste0("RSS_",arms)
  # Clear previous file
  if(nrow(results) > 100){unlink(paste0(OUTPUTS,"/results_","all.txt"))}
  # Write new file and append 
  write.table(results, file = paste0(OUTPUTS,"/results_","all.txt"),
              sep = "\t", quote = FALSE, col.names = FALSE,
              row.names = FALSE, append = TRUE)
  
  ##############
  ### RETURN ###
  
  # # Remove lines with NA
  # RSS = RSS[complete.cases(RSS),]
  
  write.table(RSS, file = paste0(wd, "/Y.txt"),
              sep = "\t", quote = FALSE, col.names = FALSE,
              row.names = FALSE)
  
  return(RSS)
  
}

# Call function
myFunction(x_parameters)


