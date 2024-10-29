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
setwd(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/round1_plots/Manuscript_Fig4"))

##########################################################################################
### Visualizes tradeoffs and reduction for hybrid vaccination and mass vaccination
##########################################################################################

data_set <- readRDS(paste0(getwd(),"/data_Fig4a.rds"))

data_set$expr = recode_factor(data_set$exp,"Moderate" ="Hybrid vaccination",
                              "ModTimed" ="Mass vaccination")

data_set <- data_set[data_set$prev_cat!="0-10" & data_set$prev_cat!="50-60" &
                       data_set$Seasons=="Short season",]

# # Plot
text_size = 18
a<-ggplot(data = data_set, aes(x = prev_cat, y = important, fill = Factors)) +
  geom_bar(position="stack", stat="identity") +
  facet_nested(outcome~expr, scales = "free_y")+
  #facet_nested(outcome~expr+Seasons, scales = "free_y")+
  scale_fill_manual(values=met.brewer("NewKingdom", 3), drop = F) +
  geom_text(aes(label = paste0(round(value,2)*100,"%")), family="Times",size = 5, hjust = 0.5,face="bold", position = position_stack(vjust = 0.5)) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),#panel.border = element_blank(),
        panel.background = element_blank(),
        #panel.grid = element_blank(),
        text=element_text(family="Times"),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0), vjust = 2, face = "bold", size=18),
        axis.title.x = element_text(margin = margin(t = 2.5, r = 0, b = 0, l = 0),vjust = 0, face = "bold", size=18),
        axis.text.x = element_text(margin = margin(t = 1), angle = 0, size=17),
        axis.text.y = element_text(margin = margin(r = 0), hjust = 0.1, size=18),
        strip.background = element_blank(),#strip.background =  element_rect(size = 0.25, colour = "black"),# strip.background = element_blank(),
        strip.text = element_text(margin = margin(t = 5), angle = 0,vjust = 1,face = "bold", size = 18),
        plot.title = element_text(hjust = 0,face = "bold", size=20),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16,face="bold"),
        legend.position = "bottom", 
        legend.direction = "horizontal")+
  labs(title = "A",  y = "Relative importance (%)", x = expression(paste(bold("Baseline "), bolditalic("Pf"), bold("PR"["2-10"]), bold("(%)"))))# +


##################
data_set <- readRDS(paste0(getwd(),"/data_Fig4b.rds"))

data_set$expr = recode_factor(data_set$exp,"Moderate" ="Hybrid vaccination",
                              "ModTimed" ="Mass vaccination")

data_set0 <- data_set[data_set$season=="Short season" &
                        data_set$out!="24m5" &
                        data_set$prev_baseline=="10-20",]
text_size = 15
b<-ggplot(data = data_set0, aes(x = x2, y = M, colour = Factors)) +
  geom_line(lwd=2)+
  geom_ribbon(aes(ymin = L, ymax = U, fill = Factors), alpha = 0.1, colour = NA)+
  facet_nested(outcome~expr, scales = "free_y")+
  #coord_cartesian(ylim = c(10,45))+
  scale_color_manual(values=met.brewer("NewKingdom", 3), drop = F) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),#panel.border = element_blank(),
        panel.background = element_blank(),
        #panel.grid = element_blank(),
        text=element_text(family="Times"),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0), vjust = 2, face = "bold", size=18),
        axis.title.x = element_text(margin = margin(t = 2.5, r = 0, b = 0, l = 0),vjust = 0, face = "bold", size=18),
        axis.text.x = element_text(margin = margin(t = 1), angle = 0, size=18),
        axis.text.y = element_text(margin = margin(r = 0), hjust = 0.1, size=18),
        strip.background = element_blank(),#strip.background =  element_rect(size = 0.25, colour = "black"),# strip.background = element_blank(),
        strip.text = element_text(margin = margin(t = 5), angle = 0,vjust = 1,face = "bold", size = 18),
        plot.title = element_text(hjust = 0,face = "bold", size=20),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16,face="bold"),
        legend.position = "bottom",
        legend.direction = "horizontal")+
  labs(title = "B",  y = paste0("Target Reduction (%)"), x="% of factor value")

ab = ggarrange(a,b, ncol=2, common.legend = F,legend="bottom")
ab

ggsave(filename = paste0(getwd(),"/","Fig 4.pdf"),
       plot = ab,
       width = 18,
       height = 10,
       dpi = 300)
