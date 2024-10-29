###################
##Plotting vaccine decay over 2 years with boosters
###################
rm(list=ls())

library(dplyr)

vaccine_decay = data.frame(matrix(nrow = 25))
vaccine_decay[,1] = seq(0:24) - 1
names(vaccine_decay) = "mos"

#91.1% [95% CI 74.5–99.7%] 

# weibull function known for RTS,S primary series
L= 7.3; k = 0.69; eff_init = 91.1; eff_init_lb = 74.5; eff_init_ub = 99.7 

# decay following first dose
vaccine_decay$eff_primary = eff_init * exp( -(vaccine_decay$mos/L)^k * log(2) )
vaccine_decay$eff_primary_lb = eff_init_lb * exp( -(vaccine_decay$mos/L)^k * log(2) )
vaccine_decay$eff_primary_ub = eff_init_ub * exp( -(vaccine_decay$mos/L)^k * log(2) )

# weibull function estimated for RTS,S after booster dose
L= 7.3; k = 0.69; eff_init_boost = 75.0

# decay following booster
vaccine_decay$mos_boost = vaccine_decay$mos - 11
vaccine_decay$eff_booster = NA
vaccine_decay[c(12:25),"eff_booster"] = eff_init_boost * exp( -(vaccine_decay[c(12:25),"mos_boost"]/L)^k * log(2) )

vaccine_decay$eff_booster[is.na(vaccine_decay$eff_booster)] = vaccine_decay$eff_primary[is.na(vaccine_decay$eff_booster)]

# plot vaccine decay
a = ggplot(data=vaccine_decay) +
  geom_line(aes(x=mos,y=eff_booster, color="with Booster"), size=1.8) +
  geom_line(aes(x=mos,y=eff_primary, color="without Booster"), size=1.4, lty = 2) +
  geom_ribbon(mapping = aes(x = mos, ymin = eff_primary_lb, ymax = eff_primary_ub),
              data = vaccine_decay, linetype=2, alpha=0.1,fill="deepskyblue2")+
  coord_cartesian(ylim = c(0.0,100), xlim=c(0,24)) +
  scale_x_continuous(breaks = seq(0,24,6))+
  scale_color_manual(values = c('without Booster' = 'deepskyblue4', 'with Booster' = 'lightskyblue4')) +
  labs(color = 'Vaccine decay \n      profiles')+
  labs(x = "Time since dose 3 (months)", title = "A", y = ("Initial efficacy (%)"))+
  theme(text=element_text(family="Times"),
        panel.border=element_rect(linetype=1,fill=NA),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 18),
        axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size = 14,face="bold"),
        legend.title =element_text(size = 14,face="bold"),
        legend.text = element_text(size=14), 
        legend.position = c(0.225,0.2))

# Seasonal profiles used
seasonality = read.table(paste0("~/M3TPP/Experiments/round1_plots/S1 Text_Fig B/Seasonality_BF_Mali.txt"), sep="\t", header = TRUE)
seasonality = gather(seasonality, mos, value, Jan:Dec, factor_key=TRUE)
seasonality$Country = recode_factor(seasonality$Seasonality,
                                    "BF" = "Burkina Faso", "Mali" = "Mali")

# plot vaccine decay
b = ggplot(data=seasonality) +
  geom_line(aes(x=mos,y=value, group=Country, color=Country), size=1.8) +
  coord_cartesian(ylim = c(0.0,40)) +
  scale_color_manual(values = c('Burkina Faso' = 'seagreen4', 'Mali' = 'seagreen2')) +
  labs(x = "Time", title = "B", 
       y = ("Proportion of annual transmission (%)"), color = 'Country')+
  theme(text=element_text(family="Times"),
        panel.border=element_rect(linetype=1,fill=NA),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 16),
        axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size = 14,face="bold"),
        legend.title =element_text(size = 14,face="bold"),
        legend.text = element_text(size=14), 
        legend.position = c(0.2,0.2))

ab = ggarrange(a,b, ncol=2, common.legend = F)
ab

ggsave(filename = paste0(getwd(),"/","Fig B.pdf"),
       plot = ab,
       width = 12,
       height = 6,
       dpi = 300)
