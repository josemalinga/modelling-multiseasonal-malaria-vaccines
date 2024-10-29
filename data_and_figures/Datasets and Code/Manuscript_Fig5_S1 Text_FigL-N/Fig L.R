rm(list=ls())
# Library
library(ggplot2)
library(dplyr)
library(hetGP)
library(MetBrewer)
library(gmodels)
library("ggsci")
library("wesanderson")
#names(wes_palettes)
library(ggh4x)
library(ggpubr)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory
setwd(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/round1_plots/Manuscript_Fig5"))

##########################################################################################
### Visualizes tradeoffs and reduction and sensitivity analysis for hybrid vaccination
### Multiseasonal impact
##########################################################################################

data_set <- readRDS(paste0(getwd(),"/data_multiseasonal_Tradeoffs.rds"))

data_set$expr = recode_factor(data_set$exp,"multi" ="Hybrid vaccination",
                              "multiBoost" ="Hybrid vaccination",
                              "multiTimed" ="Mass vaccination",
                              "multiTimedBoost" ="Mass vaccination",)

data_set$boost = recode_factor(data_set$exp,"multi" ="without Booster",
                               "multiBoost" ="with Booster",
                               "multiTimed" ="without Booster",
                               "multiTimedBoost" ="with Booster",)

data_set <- data_set[data_set$expr=="Hybrid vaccination",]

data_set$Seasons = recode_factor(data_set$Seasons,"Seasonal transmission" ="Short season",
                                 "Long season" ="Long season","Perennial transmission"="Constant transmission" )

data_set = data_set[data_set$Coverage==1 &
                      data_set$prev_cat=="10-20" &
                      data_set$Seasons== "Long season" &
                      data_set$boost=="without Booster"
                    , ]

text_size = 15
# #when we want to check the predicted reductions by prevalence and all
a<-ggplot(data=data_set, aes(x = Duration, y = InitialEfficacy, fill = Reduction)) +
  geom_tile(width = 1.5, height = 2) +
  # facet_nested(.~expr+outcome)+
  facet_nested(.~outcome)+
  theme(panel.grid = element_blank(),
        text=element_text(family="Times"),
        panel.border = element_blank(),#panel.border = element_rect(colour = "grey", fill=NA, linewidth=0.5),
        panel.background = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0), vjust = 2, face = "bold", size=18),
        axis.title.x = element_text(margin = margin(t = -5, r = 0, b = 0, l = 0),vjust = 0, face = "bold", size=18),
        axis.text.x = element_text(margin = margin(t = 1), angle = 0, colour = "black", size=15),
        axis.text.y = element_text(margin = margin(r = 0), hjust = 0, colour = "black", size=17),
        #axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0,face = "bold", size=20),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", vjust = -0.5, size = 18),
        panel.spacing = unit(0, "lines"),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14,face="bold"),
        legend.position = "right", 
        legend.direction = "vertical",
        legend.key = element_rect(fill = "white"),
        legend.key.height = unit(1.2, 'cm'), legend.key.width= unit(1.5, 'cm'))+
  coord_cartesian(xlim = c(6,18),ylim = c(50,100))+
  scale_x_continuous(limits = c(6,18), breaks=seq(6,18,3))+
  scale_fill_stepsn(colours = c("#1B4D79", "#5880B1", "#F9A24B", "#F9D48E", "#FFEDCB"),
                    breaks= seq(0,80,5),
                    limits=c(0,80),
                    show.limits=F)+
  guides(fill = guide_colourbar(title.position = "top", order = 0,title = "  Target \nreduction (%)"))+
  labs(title = "A",  y = paste0("Initial Efficacy (%)"), x="Halflife duration (months)")#+#x = expression(paste(bold("Baseline "), bolditalic("Pf"), bold("PR"["2-10"])))) +


##########################################################################################

data_set <- readRDS(paste0(getwd(),"/data_multiseasonal_Reduction.rds"))

data_set$expr = recode_factor(data_set$exp,"multi" ="Hybrid vaccination",
                              "multiBoost" ="Hybrid vaccination",
                              "multiTimed" ="Mass vaccination",
                              "multiTimedBoost" ="Mass vaccination",)

data_set$boost = recode_factor(data_set$exp,"multi" ="without Booster",
                               "multiBoost" ="with Booster",
                               "multiTimed" ="without Booster",
                               "multiTimedBoost" ="with Booster",)

data_set <- data_set[data_set$expr=="Hybrid vaccination",]

data_set$Seasons = recode_factor(data_set$Seasons,"Seasonal transmission" ="Short season",
                                 "Long season" ="Long season","Perennial transmission"="Constant transmission" )

data_set = data_set[data_set$half_cat2!="6-12mo" &
                      data_set$prev_cat!="0-10" &
                      data_set$Seasons== "Long season" &
                      data_set$boost=="without Booster"
                    , ]

data_set$value <- data_set$Reduction
data_set0 <- data_set %>% 
  group_by(grps, prev_cat,access, Seasons, outcome,expr, boost) %>% 
  mutate(across(starts_with(c("value")),
                .fns = list(
                  #mean = ~ci(.x, na.rm=TRUE)[1],
                  #se = ~ci(.x, na.rm=TRUE)[4],
                  #n = ~nobs(.x),
                  ci_l   = ~quantile(.x, probs = 0.25),
                  ci_u   = ~quantile(.x, probs = 0.75),
                  median = ~quantile(.x, probs = 0.50)
                )))

data_set1 <- data_set %>%
  group_by(grps, prev_cat,access, Seasons, outcome,expr, boost) %>%
  reframe(value_median = mean(Reduction),sd=sd(Reduction),value_ci_u=(value_median+3*sd),value_ci_l=(value_median-3*sd))

b<-ggplot(data=data_set1) +
  geom_bar(aes(y = value_median , x = prev_cat, fill=factor(grps)), position="dodge", stat="identity")+
  geom_errorbar(aes(x=prev_cat,ymax=value_ci_u, ymin=value_ci_l, color=grps), width=.2,
                position=position_dodge(.9), show.legend = FALSE)+
  # facet_nested(.~expr+outcome)+
  facet_nested(.~outcome)+
  theme(text=element_text(family="Times"),
        panel.border = element_rect(colour = "grey", fill=NA, linewidth=0.5),#panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x.bottom = element_line(color="grey"),
        axis.line.y.left = element_line(color="grey"),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0), vjust = 2, face = "bold", size=18),
        axis.title.x = element_text(margin = margin(t = -15, r = 0, b = 0, l = 0),vjust = 0, face = "bold", size=18),
        axis.text.x = element_text(margin = margin(t = 10), angle = 45, colour = "black", size=14),
        axis.text.y = element_text(margin = margin(r = 10), hjust = 0.1, colour = "black", size=17),
        strip.text = element_text(margin = margin(t = 5), angle = 0,vjust = 1,face = "bold", size = 18),
        plot.title = element_text(hjust = 0,face = "bold", size=20),
        strip.background = element_blank(),
        # strip.text = ggtext::element_markdown(face = "bold", size = 16),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16,face="bold"),
        legend.position = "right", 
        legend.direction = "vertical",
        legend.key = element_rect(fill = "white"),
        legend.key.height = unit(1.2, 'cm'), legend.key.width= unit(0.25, 'cm'))+
  coord_cartesian(ylim = c(0,80))+
  scale_fill_manual(values = c("#2F4F4F", "#798E87", "#9C964A","#CDC08C"))+#wes_palette("Moonrise2", n = 4))+
  scale_colour_manual(values = c("#2F4F4F", "#798E87", "#9C964A","#CDC08C"))+
  labs(title = "B",fill = "Tradeoffs", y = paste0("Median relative reduction (%)"), 
       x = expression(paste(bold("Baseline "), bolditalic("Pf"), bold("PR"["2-10"]))))


##########################################################################################

data_set <- readRDS(paste0(getwd(),"/data_multiseasonal_Sensitivity.rds"))

data_set$expr = recode_factor(data_set$exp,"multi" ="Hybrid vaccination",
                              "multiBoost" ="Hybrid vaccination",
                              "multiTimed" ="Mass vaccination",
                              "multiTimedBoost" ="Mass vaccination",)

data_set$boost = recode_factor(data_set$exp,"multi" ="without Booster",
                               "multiBoost" ="with Booster",
                               "multiTimed" ="without Booster",
                               "multiTimedBoost" ="with Booster",)

data_set <- data_set[data_set$expr=="Hybrid vaccination",]

data_set$Seasons = recode_factor(data_set$Seasons,"Seasonal transmission" ="Short season",
                                 "Long season" ="Long season","Perennial transmission"="Constant transmission" )

data_set = data_set[data_set$prev_cat!="0-10" & data_set$prev_cat!="50-60" &
                      data_set$Seasons== "Long season" &
                      data_set$boost=="without Booster" &
                      data_set$group!="Coverage", ]


data_set$Factors = recode_factor(data_set$Factors, "Coverage (100%)"="Coverage \n(100%)",
                                 "Efficacy (50-100%)"="Efficacy \n(75-100%)",
                                 "Halflife (6-18 months)"="Halflife \n(6-18 months)")

data_set$Factors <- droplevels(data_set$Factors)

text_size = 18
c<-ggplot(data = data_set, aes(x = prev_cat, y = important, fill = Factors)) +
  geom_bar(position="stack", stat="identity") +
  # facet_nested(.~expr+outcome)+
  facet_nested(.~outcome)+
  #scale_fill_manual(values=met.brewer("NewKingdom", 3), drop = F) + #"#e1846c", "#9eb4e0", "#e6bb9e", "#9c6849", "#735852"
  scale_fill_manual(values = c("#9eb4e0", "#e6bb9e"))+
  geom_text(aes(label = paste0(round(value,2)*100,"%")), family="Times",size = 5, hjust = 0.5,face="bold", position = position_stack(vjust = 0.5)) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "grey", fill=NA, size=0.5),#panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x.bottom = element_line(color="grey"),
        axis.line.y.left= element_line(color="grey"),
        text=element_text(family="Times"),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0), vjust = 2, face = "bold", size=18),
        axis.title.x = element_text(margin = margin(t = -15, r = 0, b = 0, l = 0),vjust = 0, face = "bold", size=18),
        axis.text.x = element_text(margin = margin(t = 10), angle = 45, size=14),
        axis.text.y = element_text(margin = margin(r = 0), hjust = 0.1, size=18),
        strip.background = element_blank(),#strip.background =  element_rect(size = 0.25, colour = "black"),# strip.background = element_blank(),
        strip.text = element_text(margin = margin(t = 5), angle = 0,vjust = 1,face = "bold", size = 18),
        plot.title = element_text(hjust = 0,face = "bold", size=20),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16,face="bold"),
        legend.position = "right", 
        legend.direction = "vertical",
        legend.key.height = unit(1.2, 'cm'), legend.key.width= unit(0.25, 'cm'))+
  labs(title = "C",  y = "Relative importance (%)", x = expression(paste(bold("Baseline "), bolditalic("Pf"), bold("PR"["2-10"]), bold("(%)"))))

abc = ggarrange(a,b,c, ncol=1, common.legend = F,legend="right")
abc

ggsave(filename = paste0(getwd(),"/","Fig L.pdf"),
       plot = abc,
       width = 15,
       height = 12,
       dpi = 300)


