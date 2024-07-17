#############################
# Sensitivity analysis for GP prevalence reduction 
#
# created 28.11.2018
# monica.golumbeanu@unibas.ch
#
# Updated September 2021
# lydia.braunack-mayer@swisstph.ch
#############################


user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
source(paste0("/scicore/home/penny/",user,"/M3TPP/analysisworkflow/3_GP_train/GP_toolbox.R"))

 
args = commandArgs(TRUE)
gp_file = args[1]
ranges_file = args[2]
results_folder = args[3]
scale = as.logical(args[4])
manual = as.logical(args[5])
print(paste0("results_folder:",results_folder))
print(paste0("scale arg:",scale))

exp_name = tools::file_path_sans_ext(basename(gp_file))
sidx_file = paste(results_folder, exp_name, "_sidx.RData", sep="")

# Load the GP and parameter ranges
gp_result_name = load(gp_file)
gp_result =  cv_result$GP_model
load(ranges_file)
param_ranges_cont
param_ranges_cont <- param_ranges_cont[1:3,]
param_ranges <- param_ranges_cont
param_ranges

# If inputs have been scaled to c(0, 1) when training the emulator, update parameter ranges to match
if (scale) {
  if (manual) {
    load(sub("_manual", "", ranges_file))
    param_ranges[, 1] <- (param_ranges[, 1] - param_ranges_cont[1:3, 1]) / (param_ranges_cont[1:3, 2] - param_ranges_cont[1:3, 1])
    param_ranges[, 2] <- (param_ranges[, 2] - param_ranges_cont[1:3, 1]) / (param_ranges_cont[1:3, 2] - param_ranges_cont[1:3, 1])
  } else {
    param_ranges[, 1] <- 0 
    param_ranges[, 2] <- 1 
  }
}

# Calculate the Sobol indices
sobol_idx_list = calc_sobol_idx(GP_model = gp_result, param_spec = param_ranges, num_points = 50000)
save(sobol_idx_list, file = sidx_file)