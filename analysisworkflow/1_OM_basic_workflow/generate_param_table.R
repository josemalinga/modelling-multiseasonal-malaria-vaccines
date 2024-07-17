########################################
# script generate_param_table.R
#
# construct a table with parameter values for each scenario of an OM experiment
#
# INPUTS:
#   exp: experiment name
#   param_ranges_cont: continous paramter ranges for parameters to be sampled
#   param_cat: categorical scenario parameters
#   noSamples: number of parameters sampled via lhs
#   noSeeds: no. seeds per OM simulation
#   chunk_size: batch size for simulation submission

# OUTPUTS:
#	- "big" paramter table and chunk_size parameter tables for OM submission and scenario creation in the GROUP and experiment folder

# created 21.01.2018
# monica.golumbeanu@unibas.ch
########################################
gen_paramtable <- function(exp, param_ranges_cont,param_cat, noSamples, noSeeds,chunk_size) {
  
  user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
  user_dir <- paste0(getwd(),"/Experiments/",exp)
  
  GROUP = "/scicore/home/penny/GROUP/M3TPP/"
  
  ranges_file = paste0(user_dir,"/param_ranges.RData")
  cat_file = paste0(user_dir,"/params_cat.RData")
  table_file = paste0(user_dir,"/param_tab.txt")
  scaffold_file= paste0(user_dir,"/OM_jobs/scaffold.xml")
  
  # sample continuous parameter values 
  
  Xcand = as.data.frame(lhs(noSamples, param_ranges_cont))
  colnames(Xcand) = rownames(param_ranges_cont)
  
  # "Coverage"  "Halflife"  "Efficacy"  "Efficacy2" "Coverage2"
  Xcand[,4] <- Xcand[,3]*0.75
  Xcand[,5] <- Xcand[,4]*0.75
  Xcand[,6] <- Xcand[,1]*0.8
  
  SEED = c(1:noSeeds)
  
  # save continous parameter ranges for the experiment and copy to GROUP folder
  save(param_ranges_cont, file = ranges_file)
  
  file.copy(ranges_file, paste0(GROUP,exp,"/param_ranges.RData"),overwrite=TRUE)
  
  save(param_cat, file = cat_file)
  file.copy(cat_file, paste0(GROUP,exp,"/param_ranges_cat.RData"),overwrite=TRUE)
  
  
  # Table with the parameter values
  param_all <- param_cat
  param_all[[length(param_cat)+1]] <- as.data.frame(Xcand)
  
  param_tab_all = Reduce(merge ,param_all)
  
  
  
  param_tab = cbind(scenarios_names= paste(exp, 1:nrow(param_tab_all), sep="_"),
                    param_tab_all)
  colnames(param_tab)[1] = "Scenario_Name"
  
  # Add seed column at the beginning
  param_tab = merge(param_tab, as.data.frame(SEED))
  
  # Write table to specified destination file
  file.remove(table_file)
  write.table(param_tab, table_file, sep = "\t", quote = FALSE, col.names = TRUE,
              row.names = FALSE)
  
  file.copy(table_file, paste0(GROUP,exp,"/param_tab.txt"),overwrite=TRUE)
  
  # copy over scaffold file 
  file.remove(paste0(GROUP,exp,"/scaffold.xml"))
  
  file.copy(scaffold_file, paste0(GROUP,exp,"/scaffold.xml"),overwrite=TRUE)
  
  
  
  chunk <- chunk_size
  n <- nrow(param_tab)
  r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
  split_param_tab <- split(param_tab,r)
  
  num_tab <- max(r)
  for(j in 0:(num_tab-1)){ 
    
    table_file = paste0(GROUP,exp,"/param_tab_",j,".txt")
    
    # Write table to specified destination file
    file.remove(table_file)
    write.table(split_param_tab[[j+1]], table_file, sep = "\t", quote = FALSE, col.names = TRUE,
                row.names = FALSE)
    
  }
  
}