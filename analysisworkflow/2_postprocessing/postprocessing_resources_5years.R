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
  dosetime1 = dosetime0
  
  # define age groups 
  age_groups <- c(0,0.25,0.5,0.75,1,1.75,2,2.75,3,3.75,4,4.75,5,6,7,8,9,10,15,18,55,100)
 
  #create new smaller age groups; 
  minIntAge=c(0.75,1.75,2,2.75,3.75,4.75)
  maxAge=c(1.75,2,2.75,3.75,4.75,5)
  age05 = seq(which(age_groups==0),which(age_groups==5)-1)
  age210 = seq(which(age_groups==2),which(age_groups==10)-1)
  age610 = seq(which(age_groups==6),which(age_groups==10)-1)
  age10 = seq(which(age_groups==0),which(age_groups==10)-1)
  
  # for multiseasonal vaccine 2 to 5 year olds
  # WE NEED TO FOLLOW A COHORT SO THIS FILE CHANGES
  
  # for all iTPP2 analysis for children age group aged 6 months to 5 years
  ageint = seq(which(age_groups==0),which(age_groups==6)-1)
  ageint
  
  # add age and year
  om_result$time<-as.numeric(om_result$time); om_result$age_group<-as.numeric(om_result$age_group) 
  
  # Remove first measurement as it includes all cases until then - and the first five months to make sure its 10y
  to_remove = which(as.numeric(om_result$time) ==1 )
  om_result = om_result[-to_remove,] 
  
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
  b_time = dosetime1 - (year_to_5day*2)
  btime = dosetime1 - (year_to_5day*3)
  
  dosetime5 = 675
  
  ftime_6before = btime + floor((year_to_5day/2))
  ftime_12before = btime + year_to_5day - 1
  ftime_18before = btime + year_to_5day + floor(year_to_5day/2) - 1
  ftime_24before = btime + year_to_5day*2 - 1
  
  ftime_6after5 = dosetime5 + floor(year_to_5day/2) 
  ftime_12after5 = dosetime5 + year_to_5day - 1
  ftime_18after5 = dosetime5 + year_to_5day + floor(year_to_5day/2) - 1
  ftime_24after5 = dosetime5 + (year_to_5day*2) - 1
  
  # calculate the incidence numbers for intervention ages and the average reductions
  incidence_int = nclin[nclin$age_group %in% ageint, -c(3,5)]
  incidence_int <- incidence_int %>% group_by(time) %>% 
    summarize(n = sum(n), inc=sum(inc))
  
  cpp_6before = sum(incidence_int[incidence_int$time>=btime & incidence_int$time<=ftime_6before,]$inc)/mean(incidence_int[incidence_int$time>=btime & incidence_int$time<=ftime_6before,]$n)
  cpp_12before = sum(incidence_int[incidence_int$time>=btime & incidence_int$time<=ftime_12before,]$inc)/mean(incidence_int[incidence_int$time>=btime & incidence_int$time<=ftime_12before,]$n)
  cpp_18before = sum(incidence_int[incidence_int$time>=btime & incidence_int$time<=ftime_18before,]$inc)/mean(incidence_int[incidence_int$time>=btime & incidence_int$time<=ftime_18before,]$n)
  cpp_24before = sum(incidence_int[incidence_int$time>=btime & incidence_int$time<=ftime_24before,]$inc)/mean(incidence_int[incidence_int$time>=btime & incidence_int$time<=ftime_24before,]$n)
  
  cpp_6after5 = sum(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_6after5,]$inc)/mean(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_6after5,]$n)
  cpp_12after5 = sum(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_12after5,]$inc)/mean(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_12after5,]$n)
  cpp_18after5 = sum(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_18after5,]$inc)/mean(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_18after5,]$n)
  cpp_24after5 = sum(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_24after5,]$inc)/mean(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_24after5,]$n)
  
  inc_red_int_6m5 = 100*(cpp_6before - cpp_6after5)/ cpp_6before
  inc_red_int_6m5 = inc_red_int_6m5*(inc_red_int_6m5 >= 0)
  inc_red_int_12m5 = 100*( cpp_12before - cpp_12after5 )/ cpp_12before
  inc_red_int_12m5 = inc_red_int_12m5*(inc_red_int_12m5 >= 0)
  inc_red_int_18m5 = 100*( cpp_18before - cpp_18after5 )/ cpp_18before
  inc_red_int_18m5 = inc_red_int_18m5*(inc_red_int_18m5 >= 0)
  inc_red_int_24m5 = 100*( cpp_24before - cpp_24after5 )/ cpp_24before
  inc_red_int_24m5 = inc_red_int_24m5*(inc_red_int_24m5 >= 0)
  
  tot_int_cases5y = sum(incidence_int[incidence_int$time>=dosetime1 & incidence_int$time<=(year_to_5day*10),]$inc)
  tot_int_pop5y = mean(incidence_int[incidence_int$time>=dosetime1 & incidence_int$time<=(year_to_5day*10),]$n)
  
  cpp_2nd6before = sum(incidence_int[incidence_int$time>=ftime_6before & incidence_int$time<=ftime_12before,]$inc)/mean(incidence_int[incidence_int$time>=ftime_6before & incidence_int$time<=ftime_12before,]$n)
  cpp_3rd6before = sum(incidence_int[incidence_int$time>=ftime_12before & incidence_int$time<=ftime_18before,]$inc)/mean(incidence_int[incidence_int$time>=ftime_12before & incidence_int$time<=ftime_18before,]$n)
  cpp_4th6before = sum(incidence_int[incidence_int$time>=ftime_18before & incidence_int$time<=ftime_24before,]$inc)/mean(incidence_int[incidence_int$time>=ftime_18before & incidence_int$time<=ftime_24before,]$n)
  
  cpp_2nd6after5 = sum(incidence_int[incidence_int$time>=ftime_6after5 & incidence_int$time<=ftime_12after5,]$inc)/mean(incidence_int[incidence_int$time>=ftime_6after5 & incidence_int$time<=ftime_12after5,]$n)
  cpp_3rd6after5 = sum(incidence_int[incidence_int$time>=ftime_12after5 & incidence_int$time<=ftime_18after5,]$inc)/mean(incidence_int[incidence_int$time>=ftime_12after5 & incidence_int$time<=ftime_18after5,]$n)
  cpp_4th6after5 = sum(incidence_int[incidence_int$time>=ftime_18after5 & incidence_int$time<=ftime_24after5,]$inc)/mean(incidence_int[incidence_int$time>=ftime_18after5 & incidence_int$time<=ftime_24after5,]$n)
  
  inc_red_int_2nd6m5 = 100*(cpp_2nd6before - cpp_2nd6after5)/ cpp_2nd6before
  inc_red_int_2nd6m5 = inc_red_int_2nd6m5*(inc_red_int_2nd6m5 >= 0)
 
  inc_red_int_3rd6m5 = 100*(cpp_3rd6before - cpp_3rd6after5)/ cpp_3rd6before
  inc_red_int_3rd6m5 = inc_red_int_3rd6m5*(inc_red_int_3rd6m5 >= 0)

  inc_red_int_4th6m5 = 100*(cpp_4th6before - cpp_4th6after5)/ cpp_4th6before
  inc_red_int_4th6m5 = inc_red_int_4th6m5*(inc_red_int_4th6m5 >= 0)

  
  # calculate the incidence numbers for 6-10 year olds and the average reductions
  incidence_610 = nclin[nclin$age_group %in% age610, -c(3,5)]
  incidence_610 <- incidence_610 %>% group_by(time) %>% 
    summarize(n = sum(n), inc=sum(inc))
  
  cpp_6before = sum(incidence_610[incidence_610$time>=btime & incidence_610$time<=ftime_6before,]$inc)/mean(incidence_610[incidence_610$time>=btime & incidence_610$time<=ftime_6before,]$n)
  cpp_12before = sum(incidence_610[incidence_610$time>=btime & incidence_610$time<=ftime_12before,]$inc)/mean(incidence_610[incidence_610$time>=btime & incidence_610$time<=ftime_12before,]$n)
  cpp_18before = sum(incidence_610[incidence_610$time>=btime & incidence_610$time<=ftime_18before,]$inc)/mean(incidence_610[incidence_610$time>=btime & incidence_610$time<=ftime_18before,]$n)
  cpp_24before = sum(incidence_610[incidence_610$time>=btime & incidence_610$time<=ftime_24before,]$inc)/mean(incidence_610[incidence_610$time>=btime & incidence_610$time<=ftime_24before,]$n)
  
  cpp_6after5 = sum(incidence_610[incidence_610$time>=dosetime5 & incidence_610$time<=ftime_6after5,]$inc)/mean(incidence_610[incidence_610$time>=dosetime5 & incidence_610$time<=ftime_6after5,]$n)
  cpp_12after5 = sum(incidence_610[incidence_610$time>=dosetime5 & incidence_610$time<=ftime_12after5,]$inc)/mean(incidence_610[incidence_610$time>=dosetime5 & incidence_610$time<=ftime_12after5,]$n)
  cpp_18after5 = sum(incidence_610[incidence_610$time>=dosetime5 & incidence_610$time<=ftime_18after5,]$inc)/mean(incidence_610[incidence_610$time>=dosetime5 & incidence_610$time<=ftime_18after5,]$n)
  cpp_24after5 = sum(incidence_610[incidence_610$time>=dosetime5 & incidence_610$time<=ftime_24after5,]$inc)/mean(incidence_610[incidence_610$time>=dosetime5 & incidence_610$time<=ftime_24after5,]$n)
  
  inc_red_610_6m5 = 100*(cpp_6before - cpp_6after5)/ cpp_6before
  inc_red_610_6m5 = inc_red_610_6m5*(inc_red_610_6m5 >= 0)
  inc_red_610_12m5 = 100*( cpp_12before - cpp_12after5 )/ cpp_12before
  inc_red_610_12m5 = inc_red_610_12m5*(inc_red_610_12m5 >= 0)
  inc_red_610_18m5 = 100*( cpp_18before - cpp_18after5 )/ cpp_18before
  inc_red_610_18m5 = inc_red_610_18m5*(inc_red_610_18m5 >= 0)
  inc_red_610_24m5 = 100*( cpp_24before - cpp_24after5 )/ cpp_24before
  inc_red_610_24m5 = inc_red_610_24m5*(inc_red_610_24m5 >= 0)
  
  tot_610_cases5y = sum(incidence_610[incidence_610$time>=dosetime1 & incidence_610$time<=(year_to_5day*10),]$inc)
  tot_610_pop5y = mean(incidence_610[incidence_610$time>=dosetime1 & incidence_610$time<=(year_to_5day*10),]$n)
  
  cpp_2nd6before = sum(incidence_610[incidence_610$time>=ftime_6before & incidence_610$time<=ftime_12before,]$inc)/mean(incidence_610[incidence_610$time>=ftime_6before & incidence_610$time<=ftime_12before,]$n)
  cpp_3rd6before = sum(incidence_610[incidence_610$time>=ftime_12before & incidence_610$time<=ftime_18before,]$inc)/mean(incidence_610[incidence_610$time>=ftime_12before & incidence_610$time<=ftime_18before,]$n)
  cpp_4th6before = sum(incidence_610[incidence_610$time>=ftime_18before & incidence_610$time<=ftime_24before,]$inc)/mean(incidence_610[incidence_610$time>=ftime_18before & incidence_610$time<=ftime_24before,]$n)
  
  cpp_2nd6after5 = sum(incidence_610[incidence_610$time>=ftime_6after5 & incidence_610$time<=ftime_12after5,]$inc)/mean(incidence_610[incidence_610$time>=ftime_6after5 & incidence_610$time<=ftime_12after5,]$n)
  cpp_3rd6after5 = sum(incidence_610[incidence_610$time>=ftime_12after5 & incidence_610$time<=ftime_18after5,]$inc)/mean(incidence_610[incidence_610$time>=ftime_12after5 & incidence_610$time<=ftime_18after5,]$n)
  cpp_4th6after5 = sum(incidence_610[incidence_610$time>=ftime_18after5 & incidence_610$time<=ftime_24after5,]$inc)/mean(incidence_610[incidence_610$time>=ftime_18after5 & incidence_610$time<=ftime_24after5,]$n)
  
  inc_red_610_2nd6m5 = 100*(cpp_2nd6before - cpp_2nd6after5)/ cpp_2nd6before
  inc_red_610_2nd6m5 = inc_red_610_2nd6m5*(inc_red_610_2nd6m5 >= 0)
  
  inc_red_610_3rd6m5 = 100*(cpp_3rd6before - cpp_3rd6after5)/ cpp_3rd6before
  inc_red_610_3rd6m5 = inc_red_610_3rd6m5*(inc_red_610_3rd6m5 >= 0)
  
  inc_red_610_4th6m5 = 100*(cpp_4th6before - cpp_4th6after5)/ cpp_4th6before
  inc_red_610_4th6m5 = inc_red_610_4th6m5*(inc_red_610_4th6m5 >= 0)
  
  ######################################
  #calculate output prevalence reduction
  ######################################   
  
  # Calculate the prevalence for all the monitored years 
  # Prevalence = total number of infected people (in age group)/ total population (in age group)
  prevalence_int = prev_all2[prev_all2$age_group %in% ageint, -c(3,5)]
  prevalence_int <- prevalence_int %>% group_by(time) %>% 
    summarise(n = sum(n), prev=mean(prev))
  
  cpp_6before = mean(prevalence_int[prevalence_int$time>=btime & prevalence_int$time<=ftime_6before,]$prev/prevalence_int[prevalence_int$time>=btime & prevalence_int$time<=ftime_6before,]$n)
  cpp_12before = mean(prevalence_int[prevalence_int$time>=btime & prevalence_int$time<=ftime_12before,]$prev/prevalence_int[prevalence_int$time>=btime & prevalence_int$time<=ftime_12before,]$n)
  cpp_18before = mean(prevalence_int[prevalence_int$time>=btime & prevalence_int$time<=ftime_18before,]$prev/prevalence_int[prevalence_int$time>=btime & prevalence_int$time<=ftime_18before,]$n)
  cpp_24before = mean(prevalence_int[prevalence_int$time>=btime & prevalence_int$time<=ftime_24before,]$prev/prevalence_int[prevalence_int$time>=btime & prevalence_int$time<=ftime_24before,]$n)
  
  cpp_6after5 = mean(prevalence_int[prevalence_int$time>=dosetime5 & prevalence_int$time<=ftime_6after5,]$prev/prevalence_int[prevalence_int$time>=dosetime5 & prevalence_int$time<=ftime_6after5,]$n)
  cpp_12after5 = mean(prevalence_int[prevalence_int$time>=dosetime5 & prevalence_int$time<=ftime_12after5,]$prev/prevalence_int[prevalence_int$time>=dosetime5 & prevalence_int$time<=ftime_12after5,]$n)
  cpp_18after5 = mean(prevalence_int[prevalence_int$time>=dosetime5 & prevalence_int$time<=ftime_18after5,]$prev/prevalence_int[prevalence_int$time>=dosetime5 & prevalence_int$time<=ftime_18after5,]$n)
  cpp_24after5 = mean(prevalence_int[prevalence_int$time>=dosetime5 & prevalence_int$time<=ftime_24after5,]$prev/prevalence_int[prevalence_int$time>=dosetime5 & prevalence_int$time<=ftime_24after5,]$n)
  
  prev_red_int_6m5 = 100*(cpp_6before - cpp_6after5)/ cpp_6before
  prev_red_int_6m5 = prev_red_int_6m5*(prev_red_int_6m5 >= 0)
  prev_red_int_12m5 = 100*( cpp_12before - cpp_12after5 )/ cpp_12before
  prev_red_int_12m5 = prev_red_int_12m5*(prev_red_int_12m5 >= 0)
  prev_red_int_18m5 = 100*( cpp_18before - cpp_18after5 )/ cpp_18before
  prev_red_int_18m5 = prev_red_int_18m5*(prev_red_int_18m5 >= 0)
  prev_red_int_24m5 = 100*( cpp_24before - cpp_24after5 )/ cpp_24before
  prev_red_int_24m5 = prev_red_int_24m5*(prev_red_int_24m5 >= 0)
  
  tot_int_infs5y = sum(prevalence_int[prevalence_int$time>=dosetime1 & prevalence_int$time<=(year_to_5day*10),]$prev)

  cpp_2nd6before = mean(prevalence_int[prevalence_int$time>=ftime_6before & prevalence_int$time<=ftime_12before,]$prev/prevalence_int[prevalence_int$time>=ftime_6before & prevalence_int$time<=ftime_12before,]$n)
  cpp_3rd6before = mean(prevalence_int[prevalence_int$time>=ftime_12before & prevalence_int$time<=ftime_18before,]$prev/prevalence_int[prevalence_int$time>=ftime_12before & prevalence_int$time<=ftime_18before,]$n)
  cpp_4th6before = mean(prevalence_int[prevalence_int$time>=ftime_18before & prevalence_int$time<=ftime_24before,]$prev/prevalence_int[prevalence_int$time>=ftime_18before & prevalence_int$time<=ftime_24before,]$n)
  
  cpp_2nd6after5 = mean(prevalence_int[prevalence_int$time>=ftime_6after5 & prevalence_int$time<=ftime_12after5,]$prev/prevalence_int[prevalence_int$time>=ftime_6after5 & prevalence_int$time<=ftime_12after5,]$n)
  cpp_3rd6after5 = mean(prevalence_int[prevalence_int$time>=ftime_12after5 & prevalence_int$time<=ftime_18after5,]$prev/prevalence_int[prevalence_int$time>=ftime_12after5 & prevalence_int$time<=ftime_18after5,]$n)
  cpp_4th6after5 = mean(prevalence_int[prevalence_int$time>=ftime_18after5 & prevalence_int$time<=ftime_24after5,]$prev/prevalence_int[prevalence_int$time>=ftime_18after5 & prevalence_int$time<=ftime_24after5,]$n)

  prev_red_int_2nd6m5 = 100*(cpp_2nd6before - cpp_2nd6after5)/ cpp_2nd6before
  prev_red_int_2nd6m5 = prev_red_int_2nd6m5*(prev_red_int_2nd6m5 >= 0)
  
  prev_red_int_3rd6m5 = 100*(cpp_3rd6before - cpp_3rd6after5)/ cpp_3rd6before
  prev_red_int_3rd6m5 = prev_red_int_3rd6m5*(prev_red_int_3rd6m5 >= 0)
  
  prev_red_int_4th6m5 = 100*(cpp_4th6before - cpp_4th6after5)/ cpp_4th6before
  prev_red_int_4th6m5 = prev_red_int_4th6m5*(prev_red_int_4th6m5 >= 0)
  
  prevalence_int_5y_before = mean(prevalence_int[prevalence_int$time>=btime & prevalence_int$time<dosetime1,]$prev/prevalence_int[prevalence_int$time>=btime & prevalence_int$time<dosetime1,]$n)
  
  # prevalence for all ages
  prevalence_all = prev_all[, -c(3,5)]
  prevalence_all <- prevalence_all %>% group_by(time) %>% 
    summarise(n = sum(n), prev=mean(prev))
  
  prevalence_all_5y_before = mean(prevalence_all[prevalence_all$time>=btime & prevalence_all$time<dosetime1,]$prev/prevalence_all[prevalence_all$time>=btime & prevalence_all$time<dosetime1,]$n)
  
  # prevalence in under 5 years
  # prevalence for all ages
  prevalence_u5 = prev_all[prev_all$age_group %in% age05, -c(3,5)]
  prevalence_u5 <- prevalence_u5 %>% group_by(time) %>% 
    summarise(n = sum(n), prev=mean(prev))
  
  prevalence_u5_5y_before = mean(prevalence_u5[prevalence_u5$time>=btime & prevalence_u5$time<dosetime1,]$prev/prevalence_u5[prevalence_u5$time>=btime & prevalence_u5$time<dosetime1,]$n)
  
  # Calculate the prevalence for all the monitored years 
  # Prevalence = total number of infected people (in age group)/ total population (in age group)
  prevalence_610 = prev_all2[prev_all2$age_group %in% age610, -c(3,5)]
  prevalence_610 <- prevalence_610 %>% group_by(time) %>% 
    summarise(n = sum(n), prev=mean(prev))
  
  cpp_6before = mean(prevalence_610[prevalence_610$time>=btime & prevalence_610$time<=ftime_6before,]$prev/prevalence_610[prevalence_610$time>=btime & prevalence_610$time<=ftime_6before,]$n)
  cpp_12before = mean(prevalence_610[prevalence_610$time>=btime & prevalence_610$time<=ftime_12before,]$prev/prevalence_610[prevalence_610$time>=btime & prevalence_610$time<=ftime_12before,]$n)
  cpp_18before = mean(prevalence_610[prevalence_610$time>=btime & prevalence_610$time<=ftime_18before,]$prev/prevalence_610[prevalence_610$time>=btime & prevalence_610$time<=ftime_18before,]$n)
  cpp_24before = mean(prevalence_610[prevalence_610$time>=btime & prevalence_610$time<=ftime_24before,]$prev/prevalence_610[prevalence_610$time>=btime & prevalence_610$time<=ftime_24before,]$n)
  
  cpp_6after5 = mean(prevalence_610[prevalence_610$time>=dosetime5 & prevalence_610$time<=ftime_6after5,]$prev/prevalence_610[prevalence_610$time>=dosetime5 & prevalence_610$time<=ftime_6after5,]$n)
  cpp_12after5 = mean(prevalence_610[prevalence_610$time>=dosetime5 & prevalence_610$time<=ftime_12after5,]$prev/prevalence_610[prevalence_610$time>=dosetime5 & prevalence_610$time<=ftime_12after5,]$n)
  cpp_18after5 = mean(prevalence_610[prevalence_610$time>=dosetime5 & prevalence_610$time<=ftime_18after5,]$prev/prevalence_610[prevalence_610$time>=dosetime5 & prevalence_610$time<=ftime_18after5,]$n)
  cpp_24after5 = mean(prevalence_610[prevalence_610$time>=dosetime5 & prevalence_610$time<=ftime_24after5,]$prev/prevalence_610[prevalence_610$time>=dosetime5 & prevalence_610$time<=ftime_24after5,]$n)
  
  prev_red_610_6m5 = 100*(cpp_6before - cpp_6after5)/ cpp_6before
  prev_red_610_6m5 = prev_red_610_6m5*(prev_red_610_6m5 >= 0)
  prev_red_610_12m5 = 100*( cpp_12before - cpp_12after5 )/ cpp_12before
  prev_red_610_12m5 = prev_red_610_12m5*(prev_red_610_12m5 >= 0)
  prev_red_610_18m5 = 100*( cpp_18before - cpp_18after5 )/ cpp_18before
  prev_red_610_18m5 = prev_red_610_18m5*(prev_red_610_18m5 >= 0)
  prev_red_610_24m5 = 100*( cpp_24before - cpp_24after5 )/ cpp_24before
  prev_red_610_24m5 = prev_red_610_24m5*(prev_red_610_24m5 >= 0)
  
  tot_610_infs5y = sum(prevalence_610[prevalence_610$time>=dosetime1 & prevalence_610$time<=(year_to_5day*10),]$prev)
  
  cpp_2nd6before = mean(prevalence_610[prevalence_610$time>=ftime_6before & prevalence_610$time<=ftime_12before,]$prev/prevalence_610[prevalence_610$time>=ftime_6before & prevalence_610$time<=ftime_12before,]$n)
  cpp_3rd6before = mean(prevalence_610[prevalence_610$time>=ftime_12before & prevalence_610$time<=ftime_18before,]$prev/prevalence_610[prevalence_610$time>=ftime_12before & prevalence_610$time<=ftime_18before,]$n)
  cpp_4th6before = mean(prevalence_610[prevalence_610$time>=ftime_18before & prevalence_610$time<=ftime_24before,]$prev/prevalence_610[prevalence_610$time>=ftime_18before & prevalence_610$time<=ftime_24before,]$n)
  
  cpp_2nd6after5 = mean(prevalence_610[prevalence_610$time>=ftime_6after5 & prevalence_610$time<=ftime_12after5,]$prev/prevalence_610[prevalence_610$time>=ftime_6after5 & prevalence_610$time<=ftime_12after5,]$n)
  cpp_3rd6after5 = mean(prevalence_610[prevalence_610$time>=ftime_12after5 & prevalence_610$time<=ftime_18after5,]$prev/prevalence_610[prevalence_610$time>=ftime_12after5 & prevalence_610$time<=ftime_18after5,]$n)
  cpp_4th6after5 = mean(prevalence_610[prevalence_610$time>=ftime_18after5 & prevalence_610$time<=ftime_24after5,]$prev/prevalence_610[prevalence_610$time>=ftime_18after5 & prevalence_610$time<=ftime_24after5,]$n)
  
  prev_red_610_2nd6m5 = 100*(cpp_2nd6before - cpp_2nd6after5)/ cpp_2nd6before
  prev_red_610_2nd6m5 = prev_red_610_2nd6m5*(prev_red_610_2nd6m5 >= 0)
  
  prev_red_610_3rd6m5 = 100*(cpp_3rd6before - cpp_3rd6after5)/ cpp_3rd6before
  prev_red_610_3rd6m5 = prev_red_610_3rd6m5*(prev_red_610_3rd6m5 >= 0)
  
  prev_red_610_4th6m5 = 100*(cpp_4th6before - cpp_4th6after5)/ cpp_4th6before
  prev_red_610_4th6m5 = prev_red_610_4th6m5*(prev_red_610_4th6m5 >= 0)
  
  prevalence_610_5y_before = mean(prevalence_610[prevalence_610$time>=btime & prevalence_610$time<dosetime1,]$prev/prevalence_610[prevalence_610$time>=btime & prevalence_610$time<dosetime1,]$n)
  
  # calculate the prevalence numbers for all ages and the average reductions
  prevalence_210 = prev_all[prev_all$age_group %in% age210, -c(3,5)]
  prevalence_210  <- prevalence_210  %>% group_by(time) %>% 
    summarise(n = sum(n), prev=sum(prev))

  tot_210_infs5y = sum(prevalence_210[prevalence_210$time>=dosetime1 & prevalence_210$time<=(year_to_5day*10),]$prev)
  tot_210_pop5y = mean(prevalence_210[prevalence_210$time>=dosetime1 & prevalence_210$time<=(year_to_5day*10),]$n)

    # prevalence at baseline
  prevalence_210_before = mean(prevalence_210[prevalence_210$time>=b_time & prevalence_210$time<dosetime1,]$prev/prevalence_210[prevalence_210$time>=b_time & prevalence_210$time<dosetime1,]$n)
  
  prevalence_210_5y_before = mean(prevalence_210[prevalence_210$time>=btime & prevalence_210$time<dosetime1,]$prev/prevalence_210[prevalence_210$time>=btime & prevalence_210$time<dosetime1,]$n)
  
  #####################################
  # calculate severe outcomes
  #####################################
  
  # calculate the severe numbers for all ages and the average reductions
  severe_int = nsevere[nsevere$age_group %in% ageint, -c(3,5)]
  severe_int <- severe_int %>% group_by(time) %>% 
    summarise(n = sum(n), cases=sum(cases))
  
  cpp_6before = sum(severe_int[severe_int$time>=btime & severe_int$time<=ftime_6before,]$cases)/mean(severe_int[severe_int$time>=btime & severe_int$time<=ftime_6before,]$n)
  cpp_12before = sum(severe_int[severe_int$time>=btime & severe_int$time<=ftime_12before,]$cases)/mean(severe_int[severe_int$time>=btime & severe_int$time<=ftime_12before,]$n)
  cpp_18before = sum(severe_int[severe_int$time>=btime & severe_int$time<=ftime_18before,]$cases)/mean(severe_int[severe_int$time>=btime & severe_int$time<=ftime_18before,]$n)
  cpp_24before = sum(severe_int[severe_int$time>=btime & severe_int$time<=ftime_24before,]$cases)/mean(severe_int[severe_int$time>=btime & severe_int$time<=ftime_24before,]$n)
  
  cpp_6after5 = sum(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_6after5,]$cases)/mean(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_6after5,]$n)
  cpp_12after5 = sum(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_12after5,]$cases)/mean(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_12after5,]$n)
  cpp_18after5 = sum(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_18after5,]$cases)/mean(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_18after5,]$n)
  cpp_24after5 = sum(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_24after5,]$cases)/mean(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_24after5,]$n)
  
  sev_red_int_6m5 = 100*(cpp_6before - cpp_6after5)/ cpp_6before
  sev_red_int_6m5 = sev_red_int_6m5*(sev_red_int_6m5 >= 0)
  sev_red_int_12m5 = 100*( cpp_12before - cpp_12after5 )/ cpp_12before
  sev_red_int_12m5 = sev_red_int_12m5*(sev_red_int_12m5 >= 0)
  sev_red_int_18m5 = 100*( cpp_18before - cpp_18after5 )/ cpp_18before
  sev_red_int_18m5 = sev_red_int_18m5*(sev_red_int_18m5 >= 0)
  sev_red_int_24m5 = 100*( cpp_24before - cpp_24after5 )/ cpp_24before
  sev_red_int_24m5 = sev_red_int_24m5*(sev_red_int_24m5 >= 0)
  
  tot_int_sev5y = sum(severe_int[severe_int$time>=dosetime1 & severe_int$time<=(year_to_5day*10),]$cases)
  
  cpp_2nd6before = sum(severe_int[severe_int$time>=ftime_6before & severe_int$time<=ftime_12before,]$cases)/mean(severe_int[severe_int$time>=ftime_6before & severe_int$time<=ftime_12before,]$n)
  cpp_3rd6before = sum(severe_int[severe_int$time>=ftime_12before & severe_int$time<=ftime_18before,]$cases)/mean(severe_int[severe_int$time>=ftime_12before & severe_int$time<=ftime_18before,]$n)
  cpp_4th6before = sum(severe_int[severe_int$time>=ftime_18before & severe_int$time<=ftime_24before,]$cases)/mean(severe_int[severe_int$time>=ftime_18before & severe_int$time<=ftime_24before,]$n)
  
  cpp_2nd6after5 = sum(severe_int[severe_int$time>=ftime_6after5 & severe_int$time<=ftime_12after5,]$cases)/mean(severe_int[severe_int$time>=ftime_6after5 & severe_int$time<=ftime_12after5,]$n)
  cpp_3rd6after5 = sum(severe_int[severe_int$time>=ftime_12after5 & severe_int$time<=ftime_18after5,]$cases)/mean(severe_int[severe_int$time>=ftime_12after5 & severe_int$time<=ftime_18after5,]$n)
  cpp_4th6after5 = sum(severe_int[severe_int$time>=ftime_18after5 & severe_int$time<=ftime_24after5,]$cases)/mean(severe_int[severe_int$time>=ftime_18after5 & severe_int$time<=ftime_24after5,]$n)
  
  sev_red_int_2nd6m5 = 100*(cpp_2nd6before - cpp_2nd6after5)/ cpp_2nd6before
  sev_red_int_2nd6m5 = sev_red_int_2nd6m5*(sev_red_int_2nd6m5 >= 0)
  
  sev_red_int_3rd6m5 = 100*(cpp_3rd6before - cpp_3rd6after5)/ cpp_3rd6before
  sev_red_int_3rd6m5 = sev_red_int_3rd6m5*(sev_red_int_3rd6m5 >= 0)
  
  sev_red_int_4th6m5 = 100*(cpp_4th6before - cpp_4th6after5)/ cpp_4th6before
  sev_red_int_4th6m5 = sev_red_int_4th6m5*(sev_red_int_4th6m5 >= 0)
  
  # calculate the severe numbers for all ages and the average reductions
  severe_610 = nsevere[nsevere$age_group %in% age610, -c(3,5)]
  severe_610 <- severe_610 %>% group_by(time) %>% 
    summarise(n = sum(n), cases=sum(cases))
  
  cpp_6before = sum(severe_610[severe_610$time>=btime & severe_610$time<=ftime_6before,]$cases)/mean(severe_610[severe_610$time>=btime & severe_610$time<=ftime_6before,]$n)
  cpp_12before = sum(severe_610[severe_610$time>=btime & severe_610$time<=ftime_12before,]$cases)/mean(severe_610[severe_610$time>=btime & severe_610$time<=ftime_12before,]$n)
  cpp_18before = sum(severe_610[severe_610$time>=btime & severe_610$time<=ftime_18before,]$cases)/mean(severe_610[severe_610$time>=btime & severe_610$time<=ftime_18before,]$n)
  cpp_24before = sum(severe_610[severe_610$time>=btime & severe_610$time<=ftime_24before,]$cases)/mean(severe_610[severe_610$time>=btime & severe_610$time<=ftime_24before,]$n)
  
  cpp_6after5 = sum(severe_610[severe_610$time>=dosetime5 & severe_610$time<=ftime_6after5,]$cases)/mean(severe_610[severe_610$time>=dosetime5 & severe_610$time<=ftime_6after5,]$n)
  cpp_12after5 = sum(severe_610[severe_610$time>=dosetime5 & severe_610$time<=ftime_12after5,]$cases)/mean(severe_610[severe_610$time>=dosetime5 & severe_610$time<=ftime_12after5,]$n)
  cpp_18after5 = sum(severe_610[severe_610$time>=dosetime5 & severe_610$time<=ftime_18after5,]$cases)/mean(severe_610[severe_610$time>=dosetime5 & severe_610$time<=ftime_18after5,]$n)
  cpp_24after5 = sum(severe_610[severe_610$time>=dosetime5 & severe_610$time<=ftime_24after5,]$cases)/mean(severe_610[severe_610$time>=dosetime5 & severe_610$time<=ftime_24after5,]$n)
  
  sev_red_610_6m5 = 100*(cpp_6before - cpp_6after5)/ cpp_6before
  sev_red_610_6m5 = sev_red_610_6m5*(sev_red_610_6m5 >= 0)
  sev_red_610_12m5 = 100*( cpp_12before - cpp_12after5 )/ cpp_12before
  sev_red_610_12m5 = sev_red_610_12m5*(sev_red_610_12m5 >= 0)
  sev_red_610_18m5 = 100*( cpp_18before - cpp_18after5 )/ cpp_18before
  sev_red_610_18m5 = sev_red_610_18m5*(sev_red_610_18m5 >= 0)
  sev_red_610_24m5 = 100*( cpp_24before - cpp_24after5 )/ cpp_24before
  sev_red_610_24m5 = sev_red_610_24m5*(sev_red_610_24m5 >= 0)
  
  tot_610_sev5y = sum(severe_610[severe_610$time>=dosetime1 & severe_610$time<=(year_to_5day*10),]$cases)
  
  cpp_2nd6before = sum(severe_610[severe_610$time>=ftime_6before & severe_610$time<=ftime_12before,]$cases)/mean(severe_610[severe_610$time>=ftime_6before & severe_610$time<=ftime_12before,]$n)
  cpp_3rd6before = sum(severe_610[severe_610$time>=ftime_12before & severe_610$time<=ftime_18before,]$cases)/mean(severe_610[severe_610$time>=ftime_12before & severe_610$time<=ftime_18before,]$n)
  cpp_4th6before = sum(severe_610[severe_610$time>=ftime_18before & severe_610$time<=ftime_24before,]$cases)/mean(severe_610[severe_610$time>=ftime_18before & severe_610$time<=ftime_24before,]$n)
  
  cpp_2nd6after5 = sum(severe_610[severe_610$time>=ftime_6after5 & severe_610$time<=ftime_12after5,]$cases)/mean(severe_610[severe_610$time>=ftime_6after5 & severe_610$time<=ftime_12after5,]$n)
  cpp_3rd6after5 = sum(severe_610[severe_610$time>=ftime_12after5 & severe_610$time<=ftime_18after5,]$cases)/mean(severe_610[severe_610$time>=ftime_12after5 & severe_610$time<=ftime_18after5,]$n)
  cpp_4th6after5 = sum(severe_610[severe_610$time>=ftime_18after5 & severe_610$time<=ftime_24after5,]$cases)/mean(severe_610[severe_610$time>=ftime_18after5 & severe_610$time<=ftime_24after5,]$n)
  
  sev_red_610_2nd6m5 = 100*(cpp_2nd6before - cpp_2nd6after5)/ cpp_2nd6before
  sev_red_610_2nd6m5 = sev_red_610_2nd6m5*(sev_red_610_2nd6m5 >= 0)
  
  sev_red_610_3rd6m5 = 100*(cpp_3rd6before - cpp_3rd6after5)/ cpp_3rd6before
  sev_red_610_3rd6m5 = sev_red_610_3rd6m5*(sev_red_610_3rd6m5 >= 0)
  
  sev_red_610_4th6m5 = 100*(cpp_4th6before - cpp_4th6after5)/ cpp_4th6before
  sev_red_610_4th6m5 = sev_red_610_4th6m5*(sev_red_610_4th6m5 >= 0)
  
  #####################################
  # calculate mortality/deaths outcomes
  #####################################
  
  # calculate the incidence numbers for all ages and the average reductions
  deaths_int = ndeaths[ndeaths$age_group %in% ageint, -c(3,5)]
  deaths_int <- deaths_int %>% group_by(time) %>% 
    summarise(n = sum(n), deaths=sum(deaths))
  
  cpp_6before = sum(deaths_int[deaths_int$time>=btime & deaths_int$time<=ftime_6before,]$deaths)/mean(deaths_int[deaths_int$time>=btime & deaths_int$time<=ftime_6before,]$n)
  cpp_12before = sum(deaths_int[deaths_int$time>=btime & deaths_int$time<=ftime_12before,]$deaths)/mean(deaths_int[deaths_int$time>=btime & deaths_int$time<=ftime_12before,]$n)
  cpp_18before = sum(deaths_int[deaths_int$time>=btime & deaths_int$time<=ftime_18before,]$deaths)/mean(deaths_int[deaths_int$time>=btime & deaths_int$time<=ftime_18before,]$n)
  cpp_24before = sum(deaths_int[deaths_int$time>=btime & deaths_int$time<=ftime_24before,]$deaths)/mean(deaths_int[deaths_int$time>=btime & deaths_int$time<=ftime_24before,]$n)
  
  cpp_6after5 = sum(deaths_int[deaths_int$time>=dosetime5 & deaths_int$time<=ftime_6after5,]$deaths)/mean(deaths_int[deaths_int$time>=dosetime5 & deaths_int$time<=ftime_6after5,]$n)
  cpp_12after5 = sum(deaths_int[deaths_int$time>=dosetime5 & deaths_int$time<=ftime_12after5,]$deaths)/mean(deaths_int[deaths_int$time>=dosetime5 & deaths_int$time<=ftime_12after5,]$n)
  cpp_18after5 = sum(deaths_int[deaths_int$time>=dosetime5 & deaths_int$time<=ftime_18after5,]$deaths)/mean(deaths_int[deaths_int$time>=dosetime5 & deaths_int$time<=ftime_18after5,]$n)
  cpp_24after5 = sum(deaths_int[deaths_int$time>=dosetime5 & deaths_int$time<=ftime_24after5,]$deaths)/mean(deaths_int[deaths_int$time>=dosetime5 & deaths_int$time<=ftime_24after5,]$n)
  
  dea_red_int_6m5 = 100*(cpp_6before - cpp_6after5)/ cpp_6before
  dea_red_int_6m5 = dea_red_int_6m5*(dea_red_int_6m5 >= 0)
  dea_red_int_12m5 = 100*( cpp_12before - cpp_12after5 )/ cpp_12before
  dea_red_int_12m5 = dea_red_int_12m5*(dea_red_int_12m5 >= 0)
  dea_red_int_18m5 = 100*( cpp_18before - cpp_18after5 )/ cpp_18before
  dea_red_int_18m5 = dea_red_int_18m5*(dea_red_int_18m5 >= 0)
  dea_red_int_24m5 = 100*( cpp_24before - cpp_24after5 )/ cpp_24before
  dea_red_int_24m5 = dea_red_int_24m5*(dea_red_int_24m5 >= 0)
  
  tot_int_dea5y = sum(deaths_int[deaths_int$time>=dosetime1 & deaths_int$time<=(year_to_5day*10),]$deaths)
  
  cpp_2nd6before = sum(deaths_int[deaths_int$time>=ftime_6before & deaths_int$time<=ftime_12before,]$deaths)/mean(deaths_int[deaths_int$time>=ftime_6before & deaths_int$time<=ftime_12before,]$n)
  cpp_3rd6before = sum(deaths_int[deaths_int$time>=ftime_12before & deaths_int$time<=ftime_18before,]$deaths)/mean(deaths_int[deaths_int$time>=ftime_12before & deaths_int$time<=ftime_18before,]$n)
  cpp_4th6before = sum(deaths_int[deaths_int$time>=ftime_18before & deaths_int$time<=ftime_24before,]$deaths)/mean(deaths_int[deaths_int$time>=ftime_18before & deaths_int$time<=ftime_24before,]$n)
  
  cpp_2nd6after5 = sum(deaths_int[deaths_int$time>=ftime_6after5 & deaths_int$time<=ftime_12after5,]$deaths)/mean(deaths_int[deaths_int$time>=ftime_6after5 & deaths_int$time<=ftime_12after5,]$n)
  cpp_3rd6after5 = sum(deaths_int[deaths_int$time>=ftime_12after5 & deaths_int$time<=ftime_18after5,]$deaths)/mean(deaths_int[deaths_int$time>=ftime_12after5 & deaths_int$time<=ftime_18after5,]$n)
  cpp_4th6after5 = sum(deaths_int[deaths_int$time>=ftime_18after5 & deaths_int$time<=ftime_24after5,]$deaths)/mean(deaths_int[deaths_int$time>=ftime_18after5 & deaths_int$time<=ftime_24after5,]$n)
  
  dea_red_int_2nd6m5 = 100*(cpp_2nd6before - cpp_2nd6after5)/ cpp_2nd6before
  dea_red_int_2nd6m5 = dea_red_int_2nd6m5*(dea_red_int_2nd6m5 >= 0)
  
  dea_red_int_3rd6m5 = 100*(cpp_3rd6before - cpp_3rd6after5)/ cpp_3rd6before
  dea_red_int_3rd6m5 = dea_red_int_3rd6m5*(dea_red_int_3rd6m5 >= 0)
  
  dea_red_int_4th6m5 = 100*(cpp_4th6before - cpp_4th6after5)/ cpp_4th6before
  dea_red_int_4th6m5 = dea_red_int_4th6m5*(dea_red_int_4th6m5 >= 0)
  
  # calculate the incidence numbers for all ages and the average reductions
  deaths_610 = ndeaths[ndeaths$age_group %in% age610, -c(3,5)]
  deaths_610 <- deaths_610 %>% group_by(time) %>% 
    summarise(n = sum(n), deaths=sum(deaths))
  
  cpp_6before = sum(deaths_610[deaths_610$time>=btime & deaths_610$time<=ftime_6before,]$deaths)/mean(deaths_610[deaths_610$time>=btime & deaths_610$time<=ftime_6before,]$n)
  cpp_12before = sum(deaths_610[deaths_610$time>=btime & deaths_610$time<=ftime_12before,]$deaths)/mean(deaths_610[deaths_610$time>=btime & deaths_610$time<=ftime_12before,]$n)
  cpp_18before = sum(deaths_610[deaths_610$time>=btime & deaths_610$time<=ftime_18before,]$deaths)/mean(deaths_610[deaths_610$time>=btime & deaths_610$time<=ftime_18before,]$n)
  cpp_24before = sum(deaths_610[deaths_610$time>=btime & deaths_610$time<=ftime_24before,]$deaths)/mean(deaths_610[deaths_610$time>=btime & deaths_610$time<=ftime_24before,]$n)
  
  cpp_6after5 = sum(deaths_610[deaths_610$time>=dosetime5 & deaths_610$time<=ftime_6after5,]$deaths)/mean(deaths_610[deaths_610$time>=dosetime5 & deaths_610$time<=ftime_6after5,]$n)
  cpp_12after5 = sum(deaths_610[deaths_610$time>=dosetime5 & deaths_610$time<=ftime_12after5,]$deaths)/mean(deaths_610[deaths_610$time>=dosetime5 & deaths_610$time<=ftime_12after5,]$n)
  cpp_18after5 = sum(deaths_610[deaths_610$time>=dosetime5 & deaths_610$time<=ftime_18after5,]$deaths)/mean(deaths_610[deaths_610$time>=dosetime5 & deaths_610$time<=ftime_18after5,]$n)
  cpp_24after5 = sum(deaths_610[deaths_610$time>=dosetime5 & deaths_610$time<=ftime_24after5,]$deaths)/mean(deaths_610[deaths_610$time>=dosetime5 & deaths_610$time<=ftime_24after5,]$n)
  
  dea_red_610_6m5 = 100*(cpp_6before - cpp_6after5)/ cpp_6before
  dea_red_610_6m5 = dea_red_610_6m5*(dea_red_610_6m5 >= 0)
  dea_red_610_12m5 = 100*( cpp_12before - cpp_12after5 )/ cpp_12before
  dea_red_610_12m5 = dea_red_610_12m5*(dea_red_610_12m5 >= 0)
  dea_red_610_18m5 = 100*( cpp_18before - cpp_18after5 )/ cpp_18before
  dea_red_610_18m5 = dea_red_610_18m5*(dea_red_610_18m5 >= 0)
  dea_red_610_24m5 = 100*( cpp_24before - cpp_24after5 )/ cpp_24before
  dea_red_610_24m5 = dea_red_610_24m5*(dea_red_610_24m5 >= 0)
  
  tot_610_dea5y = sum(deaths_610[deaths_610$time>=dosetime1 & deaths_610$time<=(year_to_5day*10),]$deaths)
  
  cpp_2nd6before = sum(deaths_610[deaths_610$time>=ftime_6before & deaths_610$time<=ftime_12before,]$deaths)/mean(deaths_610[deaths_610$time>=ftime_6before & deaths_610$time<=ftime_12before,]$n)
  cpp_3rd6before = sum(deaths_610[deaths_610$time>=ftime_12before & deaths_610$time<=ftime_18before,]$deaths)/mean(deaths_610[deaths_610$time>=ftime_12before & deaths_610$time<=ftime_18before,]$n)
  cpp_4th6before = sum(deaths_610[deaths_610$time>=ftime_18before & deaths_610$time<=ftime_24before,]$deaths)/mean(deaths_610[deaths_610$time>=ftime_18before & deaths_610$time<=ftime_24before,]$n)
  
  cpp_2nd6after5 = sum(deaths_610[deaths_610$time>=ftime_6after5 & deaths_610$time<=ftime_12after5,]$deaths)/mean(deaths_610[deaths_610$time>=ftime_6after5 & deaths_610$time<=ftime_12after5,]$n)
  cpp_3rd6after5 = sum(deaths_610[deaths_610$time>=ftime_12after5 & deaths_610$time<=ftime_18after5,]$deaths)/mean(deaths_610[deaths_610$time>=ftime_12after5 & deaths_610$time<=ftime_18after5,]$n)
  cpp_4th6after5 = sum(deaths_610[deaths_610$time>=ftime_18after5 & deaths_610$time<=ftime_24after5,]$deaths)/mean(deaths_610[deaths_610$time>=ftime_18after5 & deaths_610$time<=ftime_24after5,]$n)
  
  dea_red_610_2nd6m5 = 100*(cpp_2nd6before - cpp_2nd6after5)/ cpp_2nd6before
  dea_red_610_2nd6m5 = dea_red_610_2nd6m5*(dea_red_610_2nd6m5 >= 0)
  
  dea_red_610_3rd6m5 = 100*(cpp_3rd6before - cpp_3rd6after5)/ cpp_3rd6before
  dea_red_610_3rd6m5 = dea_red_610_3rd6m5*(dea_red_610_3rd6m5 >= 0)
  
  dea_red_610_4th6m5 = 100*(cpp_4th6before - cpp_4th6after5)/ cpp_4th6before
  dea_red_610_4th6m5 = dea_red_610_4th6m5*(dea_red_610_4th6m5 >= 0)
  
    ######################################
  # return results
  ######################################
  
  if(cont==FALSE) {
    # Final row with outputs to return
    return_row = cbind.data.frame(scenario_params$Scenario_Name,scenario_params$SEED,prevalence_210_before,prevalence_210_5y_before,
                                  prevalence_int_5y_before,prevalence_610_5y_before,prevalence_all_5y_before,
                                  tot_int_pop5y,tot_int_cases5y,tot_int_infs5y,tot_int_sev5y,tot_int_dea5y,tot_210_pop5y,tot_210_infs5y,
                                  tot_610_pop5y,tot_610_cases5y,tot_610_infs5y,tot_610_sev5y,tot_610_dea5y,
                                  inc_red_int_6m5,inc_red_int_12m5,inc_red_int_18m5,inc_red_int_24m5,inc_red_int_2nd6m5,inc_red_int_3rd6m5,inc_red_int_4th6m5,
                                  sev_red_int_6m5,sev_red_int_12m5,sev_red_int_18m5,sev_red_int_24m5,sev_red_int_2nd6m5,sev_red_int_3rd6m5,sev_red_int_4th6m5,
                                  dea_red_int_6m5,dea_red_int_12m5,dea_red_int_18m5,dea_red_int_24m5,dea_red_int_2nd6m5,dea_red_int_3rd6m5,dea_red_int_4th6m5,
                                  inc_red_610_6m5,inc_red_610_12m5,inc_red_610_18m5,inc_red_610_24m5,inc_red_610_2nd6m5,inc_red_610_3rd6m5,inc_red_610_4th6m5,
                                  sev_red_610_6m5,sev_red_610_12m5,sev_red_610_18m5,sev_red_610_24m5,sev_red_610_2nd6m5,sev_red_610_3rd6m5,sev_red_610_4th6m5,
                                  dea_red_610_6m5,dea_red_610_12m5,dea_red_610_18m5,dea_red_610_24m5,dea_red_610_2nd6m5,dea_red_610_3rd6m5,dea_red_610_4th6m5,
                                  prev_red_int_6m5,prev_red_int_12m5,prev_red_int_18m5,prev_red_int_24m5,prev_red_int_2nd6m5,prev_red_int_3rd6m5,prev_red_int_4th6m5,
                                  prev_red_610_6m5,prev_red_610_12m5,prev_red_610_18m5,prev_red_610_24m5,prev_red_610_2nd6m5,prev_red_610_3rd6m5,prev_red_610_4th6m5    )
    
    
    colnames(return_row) = c("Scenario_Name","seed","prevalence_210_before","prevalence_210_5y_before",
                             "prevalence_int_5y_before","prevalence_610_5y_before","prevalence_all_5y_before",
                             "tot_int_pop5y","tot_int_cases5y","tot_int_infs5y","tot_int_sev5y","tot_int_dea5y","tot_210_pop5y","tot_210_infs5y",
                             "tot_610_pop5y","tot_610_cases5y","tot_610_infs5y","tot_610_sev5y","tot_610_dea5y",
                             "inc_red_int_6m5","inc_red_int_12m5","inc_red_int_18m5","inc_red_int_24m5","inc_red_int_2nd6m5","inc_red_int_3rd6m5","inc_red_int_4th6m5",
                             "sev_red_int_6m5","sev_red_int_12m5","sev_red_int_18m5","sev_red_int_24m5","sev_red_int_2nd6m5","sev_red_int_3rd6m5","sev_red_int_4th6m5",
                             "dea_red_int_6m5","dea_red_int_12m5","dea_red_int_18m5","dea_red_int_24m5","dea_red_int_2nd6m5","dea_red_int_3rd6m5","dea_red_int_4th6m5",
                             "inc_red_610_6m5","inc_red_610_12m5","inc_red_610_18m5","inc_red_610_24m5","inc_red_610_2nd6m5","inc_red_610_3rd6m5","inc_red_610_4th6m5",
                             "sev_red_610_6m5","sev_red_610_12m5","sev_red_610_18m5","sev_red_610_24m5","sev_red_610_2nd6m5","sev_red_610_3rd6m5","sev_red_610_4th6m5",
                             "dea_red_610_6m5","dea_red_610_12m5","dea_red_610_18m5","dea_red_610_24m5","dea_red_610_2nd6m5","dea_red_610_3rd6m5","dea_red_610_4th6m5",
                             "prev_red_int_6m5","prev_red_int_12m5","prev_red_int_18m5","prev_red_int_24m5","prev_red_int_2nd6m5","prev_red_int_3rd6m5","prev_red_int_4th6m5",
                             "prev_red_610_6m5","prev_red_610_12m5","prev_red_610_18m5","prev_red_610_24m5","prev_red_610_2nd6m5","prev_red_610_3rd6m5","prev_red_610_4th6m5"
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
    return_row = cbind.data.frame(scenario_params$Scenario_Name,scenario_params$SEED,prevalence_210_before,prevalence_210_5y_before,
                                  prevalence_int_5y_before,prevalence_610_5y_before,prevalence_all_5y_before,
                                  tot_int_pop5y,tot_int_cases5y,tot_int_infs5y,tot_int_sev5y,tot_int_dea5y,tot_210_pop5y,tot_210_infs5y,
                                  tot_610_pop5y,tot_610_cases5y,tot_610_infs5y,tot_610_sev5y,tot_610_dea5y,
                                  inc_red_int_6m5,inc_red_int_12m5,inc_red_int_18m5,inc_red_int_24m5,inc_red_int_2nd6m5,inc_red_int_3rd6m5,inc_red_int_4th6m5,
                                  sev_red_int_6m5,sev_red_int_12m5,sev_red_int_18m5,sev_red_int_24m5,sev_red_int_2nd6m5,sev_red_int_3rd6m5,sev_red_int_4th6m5,
                                  dea_red_int_6m5,dea_red_int_12m5,dea_red_int_18m5,dea_red_int_24m5,dea_red_int_2nd6m5,dea_red_int_3rd6m5,dea_red_int_4th6m5,
                                  inc_red_610_6m5,inc_red_610_12m5,inc_red_610_18m5,inc_red_610_24m5,inc_red_610_2nd6m5,inc_red_610_3rd6m5,inc_red_610_4th6m5,
                                  sev_red_610_6m5,sev_red_610_12m5,sev_red_610_18m5,sev_red_610_24m5,sev_red_610_2nd6m5,sev_red_610_3rd6m5,sev_red_610_4th6m5,
                                  dea_red_610_6m5,dea_red_610_12m5,dea_red_610_18m5,dea_red_610_24m5,dea_red_610_2nd6m5,dea_red_610_3rd6m5,dea_red_610_4th6m5,
                                  prev_red_int_6m5,prev_red_int_12m5,prev_red_int_18m5,prev_red_int_24m5,prev_red_int_2nd6m5,prev_red_int_3rd6m5,prev_red_int_4th6m5,
                                  prev_red_610_6m5,prev_red_610_12m5,prev_red_610_18m5,prev_red_610_24m5,prev_red_610_2nd6m5,prev_red_610_3rd6m5,prev_red_610_4th6m5
    )
    
    colnames(return_row) = c("Scenario_Name","seed","prevalence_210_before","prevalence_210_5y_before",
                             "prevalence_int_5y_before","prevalence_610_5y_before","prevalence_all_5y_before",
                             "tot_int_pop5y","tot_int_cases5y","tot_int_infs5y","tot_int_sev5y","tot_int_dea5y","tot_210_pop5y","tot_210_infs5y",
                             "tot_610_pop5y","tot_610_cases5y","tot_610_infs5y","tot_610_sev5y","tot_610_dea5y",
                             "inc_red_int_6m5","inc_red_int_12m5","inc_red_int_18m5","inc_red_int_24m5","inc_red_int_2nd6m5","inc_red_int_3rd6m5","inc_red_int_4th6m5",
                             "sev_red_int_6m5","sev_red_int_12m5","sev_red_int_18m5","sev_red_int_24m5","sev_red_int_2nd6m5","sev_red_int_3rd6m5","sev_red_int_4th6m5",
                             "dea_red_int_6m5","dea_red_int_12m5","dea_red_int_18m5","dea_red_int_24m5","dea_red_int_2nd6m5","dea_red_int_3rd6m5","dea_red_int_4th6m5",
                             "inc_red_610_6m5","inc_red_610_12m5","inc_red_610_18m5","inc_red_610_24m5","inc_red_610_2nd6m5","inc_red_610_3rd6m5","inc_red_610_4th6m5",
                             "sev_red_610_6m5","sev_red_610_12m5","sev_red_610_18m5","sev_red_610_24m5","sev_red_610_2nd6m5","sev_red_610_3rd6m5","sev_red_610_4th6m5",
                             "dea_red_610_6m5","dea_red_610_12m5","dea_red_610_18m5","dea_red_610_24m5","dea_red_610_2nd6m5","dea_red_610_3rd6m5","dea_red_610_4th6m5",
                             "prev_red_int_6m5","prev_red_int_12m5","prev_red_int_18m5","prev_red_int_24m5","prev_red_int_2nd6m5","prev_red_int_3rd6m5","prev_red_int_4th6m5",
                             "prev_red_610_6m5","prev_red_610_12m5","prev_red_610_18m5","prev_red_610_24m5","prev_red_610_2nd6m5","prev_red_610_3rd6m5","prev_red_610_4th6m5"
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
  
  aggregated_OM =   processed_OM_sim %>% group_by(Scenario_Name) %>% summarise_at(c(names(processed_OM_sim)[which(names(processed_OM_sim)=="prevalence_210_before"):length(names(processed_OM_sim) ) ]),median,na.rm=TRUE)
  
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
