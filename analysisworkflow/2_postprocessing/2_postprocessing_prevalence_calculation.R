#######################################################
#### Assess PfPR[2-10] by EIR and seasonal profile ####
#######################################################

##############
### HEADER ###
##############
#install.packages("devtools") 
#devtools::install_github("BlakeRMills/MetBrewer") 

# Clear environment
rm(list = ls())

# Set seed for replication
set.seed(42)

# Library
library(ggplot2)
library(dplyr)
library(hetGP)
library(foreach)
library(MetBrewer)
# MetBrewer::colorblind_palettes

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory
setwd("/scicore/home/penny/GROUP/M3TPP/")

##################
### EXPERIMENT ###
##################

# Insert experiment name here
exp = "..."

# Group and simulation folder
GROUP = "/scicore/home/penny/GROUP/M3TPP/"

# Load
# read parameter table
param_table_file<- paste0(GROUP, exp, "/param_tab.txt")
scenario_params = as.data.frame(read.table(param_table_file, sep= "\t", as.is = TRUE, header = TRUE, stringsAsFactors = FALSE))

#################
### AGGREGATE ###
#################

# Get postprocessing seed file name for GP file name
seed_names = list.files(path = paste0(GROUP,exp,"/postprocessing/"), pattern = "seed", full.names = FALSE)
# seed_name = sub('.txt', '', seed_name) 
seed_names

#### Loop load data in 1 table
pp_med_results = foreach(i=1:length(seed_names), .combine = "rbind") %do% {
  pp_table = read.table(paste0(GROUP,exp,"/postprocessing/",seed_names[i]), sep= "\t", as.is = TRUE, header = TRUE, stringsAsFactors = FALSE)
  
}

summary(pp_med_results$prevalence_210_5y_before)

names(pp_med_results)
pp_med_results$prev_cat = ifelse(pp_med_results$prevalence_210_5y_before<0.10, "0-10",
                    ifelse(pp_med_results$prevalence_210_5y_before>=0.1 & pp_med_results$prevalence_210_5y_before<0.2,"10-20",
                    ifelse(pp_med_results$prevalence_210_5y_before>=0.2 & pp_med_results$prevalence_210_5y_before<0.3,"20-30",
                    ifelse(pp_med_results$prevalence_210_5y_before>=0.3 & pp_med_results$prevalence_210_5y_before<0.4,"30-40",
                    ifelse(pp_med_results$prevalence_210_5y_before>=0.4 & pp_med_results$prevalence_210_5y_before<0.5,"40-50",
                    ifelse(pp_med_results$prevalence_210_5y_before>=0.5 & pp_med_results$prevalence_210_5y_before<0.6,"50-60",
                    ifelse(pp_med_results$prevalence_210_5y_before>=0.6 & pp_med_results$prevalence_210_5y_before<0.7,"60-70", "70-80")))))))

table(pp_med_results$prev_cat, pp_med_results$EIR)
table(pp_med_results$prev_cat,pp_med_results$Seasonality)

pp_med_results$prev_210  <- factor(pp_med_results$prev_cat , levels = c("0-10","10-20","20-30","30-40","40-50"))
table(pp_med_results$prev_210,pp_med_results$Seasonality)
table(pp_med_results$prev_210, pp_med_results$EIR)


ggplot(pp_med_results, aes(as.factor(prev_210), inc_red_int_3rd6m5, fill=factor(Seasonality))) +
  facet_wrap(.~Access)+
  geom_boxplot()

pp_med_results$prev_int  <- factor(pp_med_results$prevalence , levels = c("0-10","10-20","20-30","30-40","40-50","50-60"))
table(pp_med_results$prev_int)

names(pp_med_results)
names(scenario_params)

# create a new parameter table including the baseline prevalence
names(pp_med_results)[31] <- "SEED"
sc = merge(scenario_params, pp_med_results, by=c(1:31))
sc$EIR = NULL


sc$prev_category = as.numeric(sc$prev_210)
sc$prev_category2 = as.numeric(sc$prev_int)
table(sc$prev_210)


for (i in unique(sc$prev_210)){
  for (j in unique(sc$Seasonality)){
    for (k in unique(sc$Access)){
  sc_save = sc[sc$prev_210==i & sc$Seasonality==j & sc$Access==k,]
  write.table(sc_save, paste0(GROUP,exp,"/postprocessing/","prev_",exp,"_",j,"_",i,"_",k,".txt"), sep="\t", row.names=FALSE)
    }
  }
}

names(sc)
write.table(data.frame(sc[,c(1:34,101:103)]), paste0(GROUP,exp,"/parameter_tab.txt"), sep = "\t", row.names=FALSE)




