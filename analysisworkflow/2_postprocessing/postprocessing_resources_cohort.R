##########################
# Auxiliary functions for postprocessing OpenMalaria results
#
# created 02.10.2019
# monica.golumbeanu@unibas.ch
# modified by lydia.burgert@unibas.ch 
#
#
#
# modified by josephine.malinga@unibas.ch - December 15, 2021
##########################

# library(rapportools)
# library(survival)
# library(cmprsk)
library(zoo)
#detach(package:plyr)    
library(dplyr)

# Function which calculates prevalence reduction given an OpenMalaria
# simulation result.
calculate_outputs = function(om_result, scenario_params, dosetime0, cont) {
  colnames(om_result) = c("time", "age_group", "measure", "value")
  year_to_5day = 73
  month_to_5day = 6

  # define age groups 
  age_groups <- c(0,0.25,0.5,0.75,1,1.75,2,2.75,3,3.75,4,4.75,5,6,7,8,9,10,15,18,55,100)
  
  #create new smaller age groups; 
  age210 = seq(which(age_groups==2),which(age_groups==10)-1)
  
  # add age and year
  om_result$time<-as.numeric(om_result$time); om_result$age_group<-as.numeric(om_result$age_group) 
  
  # Remove first measurement as it includes all cases until then - and the first five months to make sure its 10y
  to_remove = which(as.numeric(om_result$time) ==1 )
  om_result = om_result[-to_remove,] 
  
  # followed up cohorts (recruited at 9 months of age across )
  # om_result = om_result[om_result$age_group>=1000,]
  
  # Add date, year and month
  om_result$time = as.numeric(om_result$time)

  # get t=files by outcome, merge and calculate outcomes
  clin<-om_result[om_result$measure==14,]
  names(clin)[4]<-"inc"
  
  severe<-om_result[om_result$measure==78,]
  names(severe)[4]<-"cases"
  
  deaths<-om_result[om_result$measure==74,]
  names(deaths)[4]<-"deaths"
  
  prev<-om_result[om_result$measure==3,]
  names(prev)[4]<-"prev"
  
  # check why it appears that prevalence is lower with measure 1
  prev2<-om_result[om_result$measure==1,]
  names(prev2)[4]<-"prev"
  
  total_pop<-om_result[om_result$measure==0,]
  names(total_pop)[4]<-"n"
  
  prev_all = data.frame(inner_join(total_pop, prev, by=c("time", "age_group")))
  
  # only include simulations as children get recruited into the cohort
  total_pop = total_pop[total_pop$n>0  & om_result$age_group>=1000,]
  
  prev_all2 = data.frame(inner_join(total_pop, prev2, by=c("time", "age_group")))
  nclin =  data.frame(inner_join(total_pop, clin, by=c("time", "age_group")))
  nsevere =  data.frame(inner_join(total_pop, severe, by=c("time", "age_group")))
  ndeaths =  data.frame(inner_join(total_pop, deaths, by=c("time", "age_group")))
  
  rm(clin, prev, severe, deaths, total_pop)
  
  # sum up population intervention age-groups over the years 
  ######################################
  #calculate output clinical cases
  ######################################   
  
  # Extract the clinical case numbers 
  # set up the before and follow up years
  #dosetime1 = dosetime0
  
  # for iTPP2 main analysis, we make analysis in th 5th year of the vaccination prrogram implementation
  btime = dosetime0 - (year_to_5day*5)
  dosetime1 = btime + (year_to_5day*3)
  
  dosetime5 = dosetime0
  
  ftime_6after5 = dosetime5 + floor(year_to_5day/2) 
  ftime_12after5 = dosetime5 + year_to_5day - 1
  ftime_18after5 = dosetime5 + year_to_5day + floor(year_to_5day/2) - 1
  ftime_24after5 = dosetime5 + (year_to_5day*2) - 1
  
  # calculate the incidence numbers for intervention ages and the average reductions
  incidence_int <- nclin %>% group_by(time) %>% 
    summarize(n = sum(n), inc=sum(inc))
  
  inc_6after5 = sum(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_6after5,]$inc)
  inc_12after5 = sum(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_12after5,]$inc)
  inc_18after5 = sum(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_18after5,]$inc)
  inc_24after5 = sum(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_24after5,]$inc)
  
  inc_2nd6after5 = sum(incidence_int[incidence_int$time>=ftime_6after5 & incidence_int$time<=ftime_12after5,]$inc)
  inc_3rd6after5 = sum(incidence_int[incidence_int$time>=ftime_12after5 & incidence_int$time<=ftime_18after5,]$inc)
  inc_4th6after5 = sum(incidence_int[incidence_int$time>=ftime_18after5 & incidence_int$time<=ftime_24after5,]$inc)
  
  pop_6after5 = mean(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_6after5,]$n)
  pop_12after5 = mean(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_12after5,]$n)
  pop_18after5 = mean(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_18after5,]$n)
  pop_24after5 = mean(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_24after5,]$n)
  
  pop_2nd6after5 = mean(incidence_int[incidence_int$time>=ftime_6after5 & incidence_int$time<=ftime_12after5,]$n)
  pop_3rd6after5 = mean(incidence_int[incidence_int$time>=ftime_12after5 & incidence_int$time<=ftime_18after5,]$n)
  pop_4th6after5 = mean(incidence_int[incidence_int$time>=ftime_18after5 & incidence_int$time<=ftime_24after5,]$n)
  
  cpp_6after5_inc = sum(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_6after5,]$inc)/mean(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_6after5,]$n)
  cpp_12after5_inc = sum(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_12after5,]$inc)/mean(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_12after5,]$n)
  cpp_18after5_inc = sum(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_18after5,]$inc)/mean(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_18after5,]$n)
  cpp_24after5_inc = sum(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_24after5,]$inc)/mean(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_24after5,]$n)
  
  cpp_2nd6after5_inc = sum(incidence_int[incidence_int$time>=ftime_6after5 & incidence_int$time<=ftime_12after5,]$inc)/mean(incidence_int[incidence_int$time>=ftime_6after5 & incidence_int$time<=ftime_12after5,]$n)
  cpp_3rd6after5_inc = sum(incidence_int[incidence_int$time>=ftime_12after5 & incidence_int$time<=ftime_18after5,]$inc)/mean(incidence_int[incidence_int$time>=ftime_12after5 & incidence_int$time<=ftime_18after5,]$n)
  cpp_4th6after5_inc = sum(incidence_int[incidence_int$time>=ftime_18after5 & incidence_int$time<=ftime_24after5,]$inc)/mean(incidence_int[incidence_int$time>=ftime_18after5 & incidence_int$time<=ftime_24after5,]$n)
  
  ######################################
  #calculate output prevalence reduction
  ######################################   
  
  # Calculate the prevalence for all the monitored years 
  # Prevalence = total number of infected people (in age group)/ total population (in age group)
  prevalence_int <- prev_all2 %>% group_by(time) %>% 
    summarise(n = sum(n), prev=mean(prev))
  
  prev_6after5 = mean(prevalence_int[prevalence_int$time>=dosetime5 & prevalence_int$time<=ftime_6after5,]$prev)
  prev_12after5 = mean(prevalence_int[prevalence_int$time>=dosetime5 & prevalence_int$time<=ftime_12after5,]$prev)
  prev_18after5 = mean(prevalence_int[prevalence_int$time>=dosetime5 & prevalence_int$time<=ftime_18after5,]$prev)
  prev_24after5 = mean(prevalence_int[prevalence_int$time>=dosetime5 & prevalence_int$time<=ftime_24after5,]$prev)
  
  prev_2nd6after5 = mean(prevalence_int[prevalence_int$time>=ftime_6after5 & prevalence_int$time<=ftime_12after5,]$prev)
  prev_3rd6after5 = mean(prevalence_int[prevalence_int$time>=ftime_12after5 & prevalence_int$time<=ftime_18after5,]$prev)
  prev_4th6after5 = mean(prevalence_int[prevalence_int$time>=ftime_18after5 & prevalence_int$time<=ftime_24after5,]$prev)

  cpp_6after5_prev = mean(prevalence_int[prevalence_int$time>=dosetime5 & prevalence_int$time<=ftime_6after5,]$prev/prevalence_int[prevalence_int$time>=dosetime5 & prevalence_int$time<=ftime_6after5,]$n)
  cpp_12after5_prev = mean(prevalence_int[prevalence_int$time>=dosetime5 & prevalence_int$time<=ftime_12after5,]$prev/prevalence_int[prevalence_int$time>=dosetime5 & prevalence_int$time<=ftime_12after5,]$n)
  cpp_18after5_prev = mean(prevalence_int[prevalence_int$time>=dosetime5 & prevalence_int$time<=ftime_18after5,]$prev/prevalence_int[prevalence_int$time>=dosetime5 & prevalence_int$time<=ftime_18after5,]$n)
  cpp_24after5_prev = mean(prevalence_int[prevalence_int$time>=dosetime5 & prevalence_int$time<=ftime_24after5,]$prev/prevalence_int[prevalence_int$time>=dosetime5 & prevalence_int$time<=ftime_24after5,]$n)
  
  cpp_2nd6after5_prev = mean(prevalence_int[prevalence_int$time>=ftime_6after5 & prevalence_int$time<=ftime_12after5,]$prev/prevalence_int[prevalence_int$time>=ftime_6after5 & prevalence_int$time<=ftime_12after5,]$n)
  cpp_3rd6after5_prev = mean(prevalence_int[prevalence_int$time>=ftime_12after5 & prevalence_int$time<=ftime_18after5,]$prev/prevalence_int[prevalence_int$time>=ftime_12after5 & prevalence_int$time<=ftime_18after5,]$n)
  cpp_4th6after5_prev = mean(prevalence_int[prevalence_int$time>=ftime_18after5 & prevalence_int$time<=ftime_24after5,]$prev/prevalence_int[prevalence_int$time>=ftime_18after5 & prevalence_int$time<=ftime_24after5,]$n)

  # calculate the prevalence numbers for all ages and the average reductions
  prevalence_210 = prev_all[prev_all$age_group %in% age210, -c(3,5)]
  prevalence_210  <- prevalence_210  %>% group_by(time) %>% 
    summarise(n = sum(n), prev=sum(prev))

  # prevalence at baseline
  prevalence_210_5y_before = mean(prevalence_210[prevalence_210$time>=btime & prevalence_210$time<dosetime1,]$prev/prevalence_210[prevalence_210$time>=btime & prevalence_210$time<dosetime1,]$n)
  
  #####################################
  # calculate severe outcomes
  #####################################
  
  # calculate the severe numbers for all ages and the average reductions
  severe_int <- nsevere %>% group_by(time) %>% 
    summarise(n = sum(n), cases=sum(cases))
  
  sev_6after5 = sum(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_6after5,]$cases)
  sev_12after5 = sum(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_12after5,]$cases)
  sev_18after5 = sum(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_18after5,]$cases)
  sev_24after5 = sum(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_24after5,]$cases)
  
  sev_2nd6after5 = sum(severe_int[severe_int$time>=ftime_6after5 & severe_int$time<=ftime_12after5,]$cases)
  sev_3rd6after5 = sum(severe_int[severe_int$time>=ftime_12after5 & severe_int$time<=ftime_18after5,]$cases)
  sev_4th6after5 = sum(severe_int[severe_int$time>=ftime_18after5 & severe_int$time<=ftime_24after5,]$cases)

  cpp_6after5_sev = sum(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_6after5,]$cases)/mean(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_6after5,]$n)
  cpp_12after5_sev = sum(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_12after5,]$cases)/mean(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_12after5,]$n)
  cpp_18after5_sev = sum(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_18after5,]$cases)/mean(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_18after5,]$n)
  cpp_24after5_sev = sum(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_24after5,]$cases)/mean(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_24after5,]$n)
  
  cpp_2nd6after5_sev = sum(severe_int[severe_int$time>=ftime_6after5 & severe_int$time<=ftime_12after5,]$cases)/mean(severe_int[severe_int$time>=ftime_6after5 & severe_int$time<=ftime_12after5,]$n)
  cpp_3rd6after5_sev = sum(severe_int[severe_int$time>=ftime_12after5 & severe_int$time<=ftime_18after5,]$cases)/mean(severe_int[severe_int$time>=ftime_12after5 & severe_int$time<=ftime_18after5,]$n)
  cpp_4th6after5_sev = sum(severe_int[severe_int$time>=ftime_18after5 & severe_int$time<=ftime_24after5,]$cases)/mean(severe_int[severe_int$time>=ftime_18after5 & severe_int$time<=ftime_24after5,]$n)
  
  #####################################
  # calculate mortality/deaths outcomes
  #####################################
  
  # calculate the incidence numbers for all ages and the average reductions
  deaths_int <- ndeaths %>% group_by(time) %>% 
    summarise(n = sum(n), deaths=sum(deaths))
  
  dea_6after5 = sum(deaths_int[deaths_int$time>=dosetime5 & deaths_int$time<=ftime_6after5,]$deaths)
  dea_12after5 = sum(deaths_int[deaths_int$time>=dosetime5 & deaths_int$time<=ftime_12after5,]$deaths)
  dea_18after5 = sum(deaths_int[deaths_int$time>=dosetime5 & deaths_int$time<=ftime_18after5,]$deaths)
  dea_24after5 = sum(deaths_int[deaths_int$time>=dosetime5 & deaths_int$time<=ftime_24after5,]$deaths)
  
  dea_2nd6after5 = sum(deaths_int[deaths_int$time>=ftime_6after5 & deaths_int$time<=ftime_12after5,]$deaths)
  dea_3rd6after5 = sum(deaths_int[deaths_int$time>=ftime_12after5 & deaths_int$time<=ftime_18after5,]$deaths)
  dea_4th6after5 = sum(deaths_int[deaths_int$time>=ftime_18after5 & deaths_int$time<=ftime_24after5,]$deaths)

  cpp_6after5_dea = sum(deaths_int[deaths_int$time>=dosetime5 & deaths_int$time<=ftime_6after5,]$deaths)/mean(deaths_int[deaths_int$time>=dosetime5 & deaths_int$time<=ftime_6after5,]$n)
  cpp_12after5_dea = sum(deaths_int[deaths_int$time>=dosetime5 & deaths_int$time<=ftime_12after5,]$deaths)/mean(deaths_int[deaths_int$time>=dosetime5 & deaths_int$time<=ftime_12after5,]$n)
  cpp_18after5_dea = sum(deaths_int[deaths_int$time>=dosetime5 & deaths_int$time<=ftime_18after5,]$deaths)/mean(deaths_int[deaths_int$time>=dosetime5 & deaths_int$time<=ftime_18after5,]$n)
  cpp_24after5_dea = sum(deaths_int[deaths_int$time>=dosetime5 & deaths_int$time<=ftime_24after5,]$deaths)/mean(deaths_int[deaths_int$time>=dosetime5 & deaths_int$time<=ftime_24after5,]$n)
  
  cpp_2nd6after5_dea = sum(deaths_int[deaths_int$time>=ftime_6after5 & deaths_int$time<=ftime_12after5,]$deaths)/mean(deaths_int[deaths_int$time>=ftime_6after5 & deaths_int$time<=ftime_12after5,]$n)
  cpp_3rd6after5_dea = sum(deaths_int[deaths_int$time>=ftime_12after5 & deaths_int$time<=ftime_18after5,]$deaths)/mean(deaths_int[deaths_int$time>=ftime_12after5 & deaths_int$time<=ftime_18after5,]$n)
  cpp_4th6after5_dea = sum(deaths_int[deaths_int$time>=ftime_18after5 & deaths_int$time<=ftime_24after5,]$deaths)/mean(deaths_int[deaths_int$time>=ftime_18after5 & deaths_int$time<=ftime_24after5,]$n)
  
    ######################################
  # return results
  ######################################
  
  if(cont==FALSE) {
    # Final row with outputs to return
    return_row = cbind.data.frame(scenario_params$Scenario_Name,scenario_params$SEED,prevalence_210_5y_before,
                                  inc_6after5,inc_12after5,inc_18after5,inc_24after5,inc_2nd6after5,inc_3rd6after5,inc_4th6after5,
                                  prev_6after5,prev_12after5,prev_18after5,prev_24after5,prev_2nd6after5,prev_3rd6after5,prev_4th6after5,
                                  sev_6after5,sev_12after5,sev_18after5,sev_24after5,sev_2nd6after5,sev_3rd6after5,sev_4th6after5,
                                  dea_6after5,dea_12after5,dea_18after5,dea_24after5,dea_2nd6after5,dea_3rd6after5,dea_4th6after5,
                                  pop_6after5,pop_12after5,pop_18after5,pop_24after5,pop_2nd6after5,pop_3rd6after5,pop_4th6after5  )
    
    
    colnames(return_row) = c("Scenario_Name","seed","prevalence_210_5y_before",
                             "inc_6after5","inc_12after5","inc_18after5","inc_24after5","inc_2nd6after5","inc_3rd6after5","inc_4th6after5",
                             "prev_6after5","prev_12after5","prev_18after5","prev_24after5","prev_2nd6after5","prev_3rd6after5",'prev_4th6after5',
                             "sev_6after5","sev_12after5","sev_18after5","sev_24after5","sev_2nd6after5","sev_3rd6after5","sev_4th6after5",
                             "dea_6after5","dea_12after5","dea_18after5","dea_24after5","dea_2nd6after5","dea_3rd6after5","dea_4th6after5",
                             "pop_6after5","pop_12after5","pop_18after5","pop_24after5","pop_2nd6after5","pop_3rd6after5","pop_4th6after5"  
    )
  return(return_row)}
  else{
    out_df= list("prevalence_210"=prevalence_210,
                 "prevalence_int"=prevalence_int,
                 "incidence_int"=incidence_int,
                 "severe_int"=severe_int,
                 "deaths_int"=deaths_int
    )
    return(out_df)
  }
  if(cont==FALSE) { 
    # Final row with outputs to return
    return_row = cbind.data.frame(scenario_params$Scenario_Name,scenario_params$SEED,prevalence_210_5y_before,
                                  inc_6after5,inc_12after5,inc_18after5,inc_24after5,inc_2nd6after5,inc_3rd6after5,inc_4th6after5,
                                  prev_6after5,prev_12after5,prev_18after5,prev_24after5,prev_2nd6after5,prev_3rd6after5,prev_4th6after5,
                                  sev_6after5,sev_12after5,sev_18after5,sev_24after5,sev_2nd6after5,sev_3rd6after5,sev_4th6after5,
                                  dea_6after5,dea_12after5,dea_18after5,dea_24after5,dea_2nd6after5,dea_3rd6after5,dea_4th6after5,
                                  pop_6after5,pop_12after5,pop_18after5,pop_24after5,pop_2nd6after5,pop_3rd6after5,pop_4th6after5
    )
    
    colnames(return_row) = c("Scenario_Name","seed","prevalence_210_5y_before",
                             "inc_6after5","inc_12after5","inc_18after5","inc_24after5","inc_2nd6after5","inc_3rd6after5","inc_4th6after5",
                             "prev_6after5","prev_12after5","prev_18after5","prev_24after5","prev_2nd6after5","prev_3rd6after5",'prev_4th6after5',
                             "sev_6after5","sev_12after5","sev_18after5","sev_24after5","sev_2nd6after5","sev_3rd6after5","sev_4th6after5",
                             "dea_6after5","dea_12after5","dea_18after5","dea_24after5","dea_2nd6after5","dea_3rd6after5","dea_4th6after5",
                             "pop_6after5","pop_12after5","pop_18after5","pop_24after5","pop_2nd6after5","pop_3rd6after5","pop_4th6after5"
    )
    return(return_row)}else{
      out_df= list("prevalence_210"=prevalence_210,
                   "prevalence_int"=prevalence_int,
                   "incidence_int"=incidence_int,
                   "severe_int"=severe_int,
                   "deaths_int"=deaths_int
      )
      return(out_df)
    }
}


# Wrapper for looping across all simulation results and gathering postprocessing results in a table
postprocess_OM = function(results_folder, param_table_file, final_table_dest, 
                          final_seed_table_dest, dosetime0) {
  
  param_table = read.table(param_table_file, sep= "\t", as.is = TRUE, header = TRUE, stringsAsFactors = FALSE)
  processed_OM_sim = NULL
  for( i in 1:nrow(param_table)) {
    
    print(i)
    # Read the OM simulation result
    OM_result_file = paste(results_folder, param_table[i,]$Scenario_Name, "_", 
                           param_table[i,]$SEED, "_out.txt", sep="")
    # Calculate the necessary outputs
    if(file.exists(OM_result_file) & file.info(OM_result_file)$size > 0) {
      
      # Read in file
      OM_result = read.table(OM_result_file, sep="\t")
      
      # Identify error to skip
      mtry <- try(calculate_outputs(OM_result, param_table[i,], dosetime0,cont=FALSE),
                  silent = TRUE)
      
      # Skip error or calculate outputs
      if (class(mtry) != "try-error") {
        scenario_row = calculate_outputs(OM_result, param_table[i,], dosetime0,cont=FALSE)
        processed_OM_sim = data.frame(rbind(processed_OM_sim, scenario_row),stringsAsFactors = FALSE)
      } else {
        message("Error")
      }
    }
  }
  
  
  # Summarize results over seeds and create final results tables
  
  aggregated_OM =   processed_OM_sim %>% group_by(Scenario_Name) %>% summarise_at(c(names(processed_OM_sim)[which(names(processed_OM_sim)=="prevalence_210_5y_before"):length(names(processed_OM_sim) ) ]),median,na.rm=TRUE)
  
  no_seed_table = param_table[,-c(which(colnames(param_table)=="SEED"))]
  no_seed_table = unique(no_seed_table)
  final_seed_table = merge(no_seed_table, processed_OM_sim, by = c("Scenario_Name"))
  final_table = merge(no_seed_table, aggregated_OM, by = c("Scenario_Name"))
  
  # Write result tables (summarized and with seeds) to files
  write.table(final_table, final_table_dest, sep="\t",
              col.names = TRUE, row.names = FALSE, quote=FALSE)
  write.table(final_seed_table, final_seed_table_dest, sep="\t",
              col.names = TRUE, row.names = FALSE, quote=FALSE)
}
