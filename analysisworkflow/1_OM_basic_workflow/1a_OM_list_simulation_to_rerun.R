##############################################
##############################################
###                                        ###
### STEP 1a: LIST OF FAILED OM SIMULATIONS ###
###                                        ###
##############################################
##############################################

### -------------------------------------------------------------------------
###
### M3TPP PROJECT:
### Creates new parameter table for failed simulations
###
### adapted 10.11.2021
### narimane.nekkab@swisstph.ch
###
### adapted 10.12.2021
### josephine.malinga@swisstph.ch
###
### -------------------------------------------------------------------------

##############################
### STEP 1; RUN ON CLUSTER ###
##############################

# !!!!! Change experiment name !!!!!

# Run this code (when inside om folder): 
# cd /scicore/home/penny/GROUP/M3TPP/.../om/
# ls >> /scicore/home/penny/GROUP/M3TPP/.../filenames_of_simulations_done.txt

##############
### HEADER ###
##############

rm(list = ls())
library(sjmisc)
library(plyr)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Experiment
exp ="..."

GROUP = "/scicore/home/penny/GROUP/M3TPP/"
############
### DATA ###
############

# Get postprocessing aggregate files 
file_names = list.files(path = paste0(GROUP,exp,"/om/"), pattern = "*.txt", full.names = FALSE)
length(file_names)
write.table(file_names, paste0("/scicore/home/penny/GROUP/M3TPP/",exp, "/sim_failed2.txt"))

# Load list of successful experiments and get names files
f <- read.table(paste0("/scicore/home/penny/GROUP/M3TPP/",exp, "/sim_failed2.txt"), fill=T, header = F, stringsAsFactors = F)
f <- data.frame(f[-c(1),-c(1)]); names(f)[1]<-"V1"
f_1 <- f[grepl("_out.txt", f$V1), "V1"]
f_2 = gsub("_out.txt","", f_1)

# Load parameter table files (add more a files if more param_tab_X)
a <- read.table(paste0("/scicore/home/penny/GROUP/M3TPP/",exp,"/param_tab.txt"),header = T, stringsAsFactors = F)

# Name of simulated OM files
simulated <- which(paste0(a$Scenario_Name,"_",a$SEED) %in% f_2)

# Select failed simulations
param_tab_new <- a[-simulated,]

# Table
table(param_tab_new$Seasonality, param_tab_new$EIR)
table(a$Seasonality, a$EIR)

table( param_tab_new$EIR)
table( a$EIR)


#param_tab_new$SEED = param_tab_new$SEED+267

# Percent failures
sim_table = data.frame(Seasonality=a$Seasonality, EIR=a$EIR, value=1) 

# Save new parameter table for new simulations
write.table(param_tab_new,paste0("/scicore/home/penny/GROUP/M3TPP/",exp,"/param_tab_resubmit.txt"), sep = "\t", quote = FALSE, col.names = TRUE,
            row.names = FALSE)


# !!! Follow step 1b !!!



