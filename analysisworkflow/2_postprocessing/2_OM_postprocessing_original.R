################################
################################
###                          ###
### STEP 2: POST-PROCESSING  ###
###                          ###
################################
################################


### -------------------------------------------------------------------------
###
### M3TPP PROJECT:
### Main script for running post-processing of OM simulations 
### to aggregate data which will be used to train GP
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

# Library
#detach(package:plyr)
library(dplyr)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory
setwd(paste0("/scicore/home/penny/",user,"/M3TPP"))

# Source function scripts
source(paste0("./analysisworkflow/2_postprocessing/genOMpostprocscripts.R"))


##################
### EXPERIMENT ###
##################

# Insert experiment name here
exp ="..."

# Modelling an implementation or trial setting?
setting = "implementation"
# setting = "trial"


##################
### PARAMETERS ###
##################

# Set trial weeks (2 values; start & end week)
if(setting == "implementation"){
  trialweeks = NA
}else{
  trialweeks= ... # i.e. c(42, 66)
}

# Specify the number of years before your intervention is implemented
# Used to calculate the pre-intervention incidence and prevalence 
# yearsbeforeInt= ... 

# Specify the number of follow-up years
# Used to calculate the post-intervention prevalence and incidence 
dosetime0 = ...


##########################
### RUN POSTPROCESSING ###
##########################

# Run
genOMpostprocscripts(exp, setting, dosetime0, trialweeks)


