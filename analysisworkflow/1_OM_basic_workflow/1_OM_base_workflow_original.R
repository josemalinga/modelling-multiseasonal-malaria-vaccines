##############################
# Main script for specifying parameter values and running OpenMalaria simulations on the cluster. 
# 
#
# created 12.02.2021
# lydia.burgert@unibas.ch 
#
# updated 12.15.2021
# josephine.malinga@unibas.ch 
#############################

# Setup
rm(list = ls())
library(tgp)
set.seed(42)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory
setwd(paste0("/scicore/home/penny/",user,"/M3TPP"))

# Source function scripts
source(paste0("./analysisworkflow/1_OM_basic_workflow/genOMsimscripts.R"))
source(paste0("./analysisworkflow/1_OM_basic_workflow/generate_param_table.R"))
source(paste0("./analysisworkflow/1_OM_basic_workflow/create_folders.R"))
source(paste0("./analysisworkflow/1_OM_basic_workflow/convert_access.R"))

# Insert experiment name here
exp ="..."

# Create folder in working directory
create_folders(exp) # <----- run 1st time then comment

######################################
# Specify the desired parameter values 
######################################

chunk_size = 90000

QOS = "1day"

# Seasonality type: monthly EIR vs. Fourier transformation (choose 1)
Seasonality =read.table(paste0("./Experiments/",exp,"/seasonality_469_months_updateFeb0623.txt"), sep="\t", header = TRUE)

#############################
# Categorical variables
#############################

# Biting patterns 
Biting_pattern <- data.frame(Biting_pattern=c("ML"),
                             indoor=c(0.8),outdoor=c(0.2))


# EIR
#Moderate access
EIR= data.frame(EIR=c(0.075,0.95,1.2,3.9,7.0,10,22))

#High access = 0.7
#EIR= data.frame(EIR=c(0.075,0.95,1.2,3.9,7.0,10,22))

# age minimum and maximum
MaxAge = data.frame(MaxAge=c(4.999),maxGroup=c(8))

# Intervention decay
vdecay <- data.frame(fundecay=c("weibull"),kdecay=c(0.69),vdecay=c("exp" ) )

# Coverage of healthcare system
initial_access = c(0.3)
initial_access
Access = data.frame(Access=round(pmax(convert_access(initial_access * 100), 0.04), digits = 4))

# Combine
param_cat = list(Seasonality=Seasonality,
                 Biting_pattern=Biting_pattern,
                 EIR=EIR,
                 MaxAge=MaxAge,
                 vdecay=vdecay,
                 Access=Access)

#############################
# Continuous variables and their ranges

# Name of the experiment and parameters
param_ranges_cont = rbind( Coverage = c(0.3, 1),
                           Halflife = c(0.5, 1.5),
                           Efficacy = c(0.3, 1),
                           Efficacy2 = c(0.375,0.75),
                           Efficacy3 = c(0.375,0.75),
                           Coverage2 = c(0.4,0.8))

#param_ranges_cont = param_ranges_cont[,-2]
###############################
# Generate the parameter tables  
###############################

# no. of continuous parameters sampled via lhs
noSamples = 800

# no. of OM seeds per sample
noSeeds= 5

# Generate
#gen_paramtable(exp, param_ranges_cont,param_cat, noSamples, noSeeds,chunk_size)
gen_paramtable(exp, param_ranges_cont,param_cat, noSamples, noSeeds,chunk_size)

########################################################
# Generate submission scripts and run the OM simulations
########################################################\

# Number of outputs to get in OM folder
2*nrow(Seasonality)*nrow(Biting_pattern)*nrow(EIR)*nrow(vdecay)*nrow(MaxAge)*nrow(Access)*noSamples*noSeeds

# Run
genOMsimscripts(exp, chunk_size)



