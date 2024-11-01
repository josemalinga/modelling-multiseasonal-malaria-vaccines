##################################################
##################################################
###                                            ###
### STEP 6: GP-BASED GRID SEARCH OPTIMIZATION  ###
###                                            ###
##################################################
##################################################


### -------------------------------------------------------------------------
###
### M3TPP PROJECT:
### Additional script for running optimization procedure
### of key performance characteristics from a grid search method
### using a pre-train GP emulator
### 
### Original script:
### Created 12.02.2021
### narimane.nekkab@swisstph.ch
### lydia.braunack-mayer@swisstph.ch
### josephine.malinga@swisstph.ch
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
library(ggplot2)
library(dplyr)
library(hetGP)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory
setwd(paste0("/scicore/home/penny/",user,"/M3TPP"))


##################
### EXPERIMENT ###
##################

# Insert experiment name here
exp = "iTPP2_HybridBS_0.965_10years_ModTimed"

# Group and simulation folder
GROUP = "/scicore/home/penny/GROUP/M3TPP/"
SIM_FOLDER = paste0(GROUP,exp,"/")

##################
### PARAMETERS ###
##################

# Load parameter ranges
ranges_file = paste0(GROUP,exp,"/param_ranges.RData")
load(ranges_file)
param_ranges_cont

param_ranges_cont <- param_ranges_cont[1:3,]
param_ranges_cont

# Get postprocessing seed file name for GP file name
seed_name = list.files(path = paste0(GROUP,exp,"/postprocessing/"), pattern = "prev", full.names = FALSE)
seed_name = sub('.txt', '', seed_name) 
seed_name

# Variable to optimize
optim_param = "Halflife"

# GP using scaled values
scale = T

# Create grid of parameter values to search
ngrid <- c(8, 100, 8)
# ngrid = c("Coverage" = ..., "Halflife" = ..., "Efficacy" = ...)

# Set target range size (1 is per 1% jumps, 10 by 10% etc.)
target_range_size = 5


######################
### OUTPUT FOLDERS ###
######################

# Create output main folder and subfolder per predictor
opt_folder = paste0(SIM_FOLDER,"gp/GP_grid_optimization_red/")
if(!dir.exists(opt_folder)){
  dir.create(opt_folder)
}

# # Create error folder
# ERROR_FOLDER = paste0(GROUP,exp,"/gp/GP_grid_optimization/err/")
# if(!dir.exists(ERROR_FOLDER)){
#   dir.create(ERROR_FOLDER)
# }
# LBM: Is there anything we might want to consider writing to an error file? E.g. printing which predicted outcome is being optimised?


#######################################################
###      GRID SEARCH FUNCTION FOR OPTIMIZATION      ###
#######################################################

# Create grid search optimization function
GP_grid_search_predictions <- function(predicted, gp_file, optim_param, scale, n_grid, target_range_size, opt_folder, j){
  
  ###############
  ### SCALING ###
  ###############
  
  if(scale == TRUE){
    # Scale parameter ranges to c(0, 1)
    D <- 3
    scale_params <- t(param_ranges_cont)
    for (i in 1:D) {
      scale_params[, i] <- (scale_params[, i] - scale_params[1, i]) / (scale_params[2, i] - scale_params[1, i])
    }
  }else{
    scale_params <- t(param_ranges_cont)
  }
  
  ###################
  ### GRID SEARCH ###
  ###################
  
  # Load GP model
  gp_result_name <- load(gp_file)
  gp_result <- cv_result$GP_model
  
  # July 2023 changed results to fix one value for coverage et al

  # Create scenarios (default for 3 variables; change for specific experiments when necessary)
  scenarios <- expand.grid(seq(scale_params[1, 1], scale_params[2, 1], length.out = ngrid[1]),
                           seq(scale_params[1, 2], scale_params[2, 2], length.out = ngrid[2]),
                           seq(scale_params[1, 3], scale_params[2, 3], length.out = ngrid[3]))
  names(scenarios) <- rownames(param_ranges_cont)
  
  ####################
  ### OPTIMIZATION ###
  ####################
  
  # Make predictions using emulator
  preds <- predict(x = as.matrix(scenarios), object = gp_result)
  scenarios$mean <- preds$mean
  scenarios$sd2 <- preds$sd2
  scenarios$nugs <- preds$nugs
  
  # Covert parameter values back to original scale
  if(scale == TRUE){
    for (i in rownames(param_ranges_cont)) {
      scenarios[, i] <- scenarios[, i] * (param_ranges_cont[i, 2] - param_ranges_cont[i, 1]) + param_ranges_cont[i, 1]
    }
  }
  
  #scenarios <- scenarios[c(1601:2400, 4801:5600),]
  #scenarios <- scenarios[c(seq(3,nrow(scenarios),3)),] # CHECK AND CHANGE EVERYTIME
  
  #scenarios1 <- scenarios[(scenarios$Coverage==0.5),]
  #scenarios <- rbind(scenarios[c(seq(6,nrow(scenarios),8)),], scenarios1)
  
  #unique(scenarios$Efficacy)
  #unique(scenarios$Coverage)
  
  # Identify minimum parameter value required to reach target reduction
  scenarios$target_range <- floor(scenarios$mean/target_range_size)*target_range_size
  
  # Optimization
  optimization_results <- scenarios %>%
    group_by_at(c(rownames(param_ranges_cont)[!(rownames(param_ranges_cont) == optim_param)], "target_range")) %>%
    summarise_at(optim_param, min) %>% 
    ungroup()
  
  # Rename optimized variable
  colnames(optimization_results)[colnames(optimization_results) == optim_param] <- paste0("optimal_",optim_param)
  head(optimization_results)
  
  # Save tables
  write.table(scenarios, file=paste0(opt_folder,predicted,"/",j,"_",predicted,"_GP_grid_scenario_predictions.txt")) 
  write.table(optimization_results, file=paste0(opt_folder,predicted,"/",j,"_",predicted,"_GP_grid_optimization_",optim_param,".txt")) 
  
  return(list(scenarios, optimization_results))
}


#############################################################################
###      Performing optimization with grid search for all predictors      ###
#############################################################################

# All predictors
pred_list = c(
  # 5 year impact
  "inc_red_int_6m5", "inc_red_int_12m5", "inc_red_int_24m5",
  "prev_red_int_6m5","prev_red_int_12m5","prev_red_int_24m5",
  "sev_red_int_6m5","sev_red_int_12m5","sev_red_int_24m5",
  "dea_red_int_6m5","dea_red_int_12m5","dea_red_int_24m5"
  
  #multiseasonal impact
  #"inc_red_int_6m5","inc_red_int_2nd6m5","inc_red_int_3rd6m5","inc_red_int_4th6m5"
)

# Loop
for(i in pred_list){
  
  print(i)
  # All seeds for each predictor
  
  for (j in seed_name){
    
    print(j)
    
    # Update GP file name
    gp_file = paste0("/scicore/home/penny/GROUP/M3TPP/",exp,"/gp/trained/",i,"/",j,"_", i,"_cv.RData")
    
    # Create folder
    if(!dir.create(paste0(opt_folder,i,"/"))){
      dir.create(paste0(opt_folder,i,"/"))
    }
    
    # Run (check it works)
    results = GP_grid_search_predictions(i, gp_file, optim_param, scale, n_grid, target_range_size, opt_folder, j)
    scenario_predictions = results[[1]]
    optimization_results = results[[2]]
    #rm(results)
}
}
