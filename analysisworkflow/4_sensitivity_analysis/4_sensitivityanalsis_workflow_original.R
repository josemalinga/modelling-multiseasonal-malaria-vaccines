#####################################
#####################################
###                               ###
### STEP 4: SENSITIVITY ANALYSIS  ###
###                               ###
#####################################
#####################################


### -------------------------------------------------------------------------
###
### M3TPP PROJECT:
### Main script for running sensitivity analyses
### 
### Original script:
### Created 12.02.2021
### lydia.burgert@unibas.ch 
###
### Adapted script:
### Saved 01.09.2021
### narimane.nekkab@unibas.ch
###
### R version 3.6.0
###
### -------------------------------------------------------------------------

##############
### HEADER ###
##############

# Clear environment
rm(list = ls())

# Set seed for replication
set.seed(42)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory
setwd(paste0("/scicore/home/penny/",user,"/M3TPP"))

# Source function scripts
source(paste0("./analysisworkflow/4_sensitivity_analysis/gensensanalysisscripts.R"))


##################
### EXPERIMENT ###
##################

# Insert experiment name here
exp ="..."


##################
### PARAMETERS ###
##################

# Specify parameter to be predicted
# Make sure GP exists or will get error

# All current variables for loop --> may be updated in future
pred_list = c(
  # 5 year impact
  "inc_red_int_6m5", "inc_red_int_12m5", "inc_red_int_24m5",
  "prev_red_int_6m5","prev_red_int_12m5","prev_red_int_24m5",
  "sev_red_int_6m5","sev_red_int_12m5","sev_red_int_24m5",
  "dea_red_int_6m5","dea_red_int_12m5","dea_red_int_24m5"
  
  #multiseasonal impact
  #"inc_red_int_6m5","inc_red_int_2nd6m5","inc_red_int_3rd6m5","inc_red_int_4th6m5"
)

# Specify whether parameter ranges used for the sensitivity analysis should be specified manually, e.g. for clinical translation
manual = TRUE

# If manual = TRUE, specify ranges to use
param_ranges_cont = rbind(parameter1 = c(..., ...),
                          parameter2 = c(..., ...),
                          parameter3 = c(..., ...))
ranges_file = paste0("/scicore/home/penny/GROUP/M3TPP/", exp, "/param_ranges_manual.RData")
save(param_ranges_cont, file = ranges_file)

###############
### SCALING ###
###############

# Specify whether, when training the GP, input parameters were scaled to c(0, 1) (TRUE) or used on their original scale (FALSE)
scale = TRUE

####################
### RUN ANALYSIS ###
####################

# loop run for each
for(i in pred_list){
  predicted = i
  print(i)
  gensensanalysisscripts(exp, predicted, scale, manual)
}






