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
setwd(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/round1_plots/Manuscript_Fig3"))

##########################################################################################
### Visualizes tradeoffs and reduction for hybrid vaccination and mass vaccination
##########################################################################################

data_set <- readRDS(paste0(getwd(),"/data_Fig3a.rds"))

data_set = data_set[data_set$Coverage==0.5 &
                      data_set$prev_cat=="10-20" &
                      data_set$Seasons=="Seasonal transmission", ]

data_set$Seasons = recode_factor(data_set$Seasons,"Seasonal transmission" ="Short season",
                                 "Long season" ="Long season", "Perennial transmission"="Constant transmission")

data_set$expr = recode_factor(data_set$exp,"Moderate" ="Hybrid vaccination",
                                 "ModTimed" ="Mass vaccination")

range(data_set$Reduction)
text_size = 15

# #when we want to check the predicted reductions by prevalence and all
a<-ggplot(data=data_set, aes(x = Duration, y = InitialEfficacy, fill = Reduction)) +
  geom_tile(width = 2, height = 2) +
  facet_nested(outcome~expr) +
  theme(panel.grid = element_blank(),
        text=element_text(family="Times"),
        panel.border = element_blank(),#panel.border = element_rect(colour = "grey", fill=NA, linewidth=0.5),
        panel.background = element_blank(),
        #panel.spacing.x=unit(1, "lines"),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0), vjust = -1.2, face = "bold", size=18),
        axis.title.x = element_text(margin = margin(t = 2.5, r = 0, b = 0, l = 0),vjust = 0, face = "bold", size=18),
        axis.text.x = element_text(margin = margin(t = 1), angle = 0, colour = "black", size=17),
        axis.text.y = element_text(margin = margin(r = 0), hjust = 1.5, colour = "black", size=17),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0,face = "bold", size=20),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", vjust = -0.2, size = 18),
        panel.spacing = unit(0, "lines"),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16,face="bold"),
        legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key = element_rect(fill = "white"),
        legend.key.height = unit(0.5, 'cm'), legend.key.width= unit(2, 'cm'))+
  coord_cartesian(xlim = c(6,18),ylim = c(50,100))+
  scale_x_continuous(limits = c(6,18), breaks=seq(6,18,3))+
  scale_fill_stepsn(colours = c("#1B4D79", "#5880B1", "#F9A24B", "#F9D48E", "#FFEDCB"),
                    breaks= seq(0,60,5),
                    limits=c(0,60),
                    show.limits=F)+
  guides(fill = guide_colourbar(title.position = "left", order = 0,title = "  Target \nreduction (%)"))+
  labs(title = "A",  y = paste0("Initial Efficacy (%)"), x="Halflife duration (months)")

##################
data_set <- readRDS(paste0(getwd(),"/data_Fig3b.rds"))

data_set = data_set[data_set$half_cat2=="6-12mo" &
                      data_set$prev_cat!="0-10" &
                      data_set$Seasons== "Seasonal transmission", ]

data_set$Seasons = recode_factor(data_set$Seasons,"Seasonal transmission" ="Short season",
                                 "Long season" ="Long season", "Perennial transmission"="Constant transmission")

data_set$expr = recode_factor(data_set$exp,"Moderate" ="Hybrid vaccination",
                              "ModTimed" ="Mass vaccination")

data_set0 <- data_set %>% 
  group_by(grps, prev_cat,access, Seasons, outcome, expr) %>% 
  dplyr::mutate(across(starts_with(c("value")),
                       .fns = list(
                         #mean = ~ci(.x, na.rm=TRUE)[1],
                         #se = ~ci(.x, na.rm=TRUE)[4],
                         #n = ~nobs(.x),
                         ci_l   = ~quantile(.x, probs = 0.25),
                         ci_u   = ~quantile(.x, probs = 0.75),
                         median = ~quantile(.x, probs = 0.50)
                       )))

data_set2 <- data_set %>%
  group_by(grps, prev_cat,access, Seasons, outcome, expr) %>%
  reframe(value_median = median(Reduction),sd=sd(Reduction),value_ci_u=(value_median+1*sd),value_ci_l=(value_median-1*sd))

# #when we want to check the predicted reductions by prevalence and all
b<-ggplot(data=data_set2) +
  geom_bar(aes(y = value_median , x = prev_cat, fill=factor(grps)), position="dodge", stat="identity")+
  geom_errorbar(aes(x=prev_cat,ymax=value_ci_u, ymin=value_ci_l, color=grps), width=.2,
                position=position_dodge(.9), show.legend = FALSE)+
  facet_nested(outcome~expr)+
  theme(text=element_text(family="Times"),
        panel.border = element_rect(colour = "grey", fill=NA, linewidth=0.5),#panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x.bottom = element_line(color="grey"),
        axis.line.y.left = element_line(color="grey"),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0), vjust = 2, face = "bold", size=18),
        axis.title.x = element_text(margin = margin(t = 2.5, r = 0, b = 0, l = 0),vjust = 0, face = "bold", size=18),
        axis.text.x = element_text(margin = margin(t = 1), angle = 0, colour = "black", size=16),
        axis.text.y = element_text(margin = margin(r = 0), hjust = 0.1, colour = "black", size=17),
        strip.text = element_text(margin = margin(t = 5), angle = 0,vjust = 1,face = "bold", size = 18),
        plot.title = element_text(hjust = 0,face = "bold", size=20),
        strip.background = element_blank(),
        # strip.text = ggtext::element_markdown(face = "bold", size = 16),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16,face="bold"),
        legend.position = "bottom", 
        legend.direction = "horizontal")+
  coord_cartesian(ylim = c(0,50))+
  scale_fill_manual(values = c("#2F4F4F", "#798E87", "#9C964A","#CDC08C"))+#wes_palette("Moonrise2", n = 4))+
  scale_colour_manual(values = c("#2F4F4F", "#798E87", "#9C964A","#CDC08C"))+
  labs(title = "B",fill = "Tradeoffs", y = paste0("Median relative reduction (%)"), 
       x = expression(paste(bold("Baseline "), bolditalic("Pf"), bold("PR"["2-10"]))))

ab = ggarrange(a,b, ncol=2, common.legend = F,legend="bottom")
ab

ggsave(filename = paste0(getwd(),"/","Fig H.pdf"),
       plot = ab,
       width = 18,
       height = 10,
       dpi = 300)