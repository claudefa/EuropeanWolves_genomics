# Correlations het, load, froh modern wolves
library(ggplot2)
library(scales)
library(statebins)
library("ggpubr")
library(tidyr)
library(readxl)
library(ggrepel)
library(ggConvexHull)
library(bedtoolsr)
library(stringr)

setwd("~/Documents/OneDrive - University of Copenhagen/ResQ_MSCA/")

#Metadata
metadata_mw <- read_excel("Files/Metadata_available/Metadata_modernpaper.xlsx", sheet = 1)
samples <- metadata_mw$SampleID

groups <- read_excel("Files/Metadata_available/Metadata_modernpaper.xlsx", sheet = 2)

order_popid <- c("EASIA","VOUR","SWASIA","KAR",
                 "EFIN","WFIN", "SCAN", "NBAL",
                 "SBAL","CAR" ,"EEL","CEL",
                 "DINBAL","NWIB","ITA")

order_pop <- c("Karelian","Eastern Finland",
               "Western Finland", "Scandinavian", "Northern Baltic",
               "Southern Baltic","Carpathian" ,"Eastern European Lowlands","Central European Lowlands",
               "Dinaric Balkan","North Western Iberian","Italian",
               "Eastern Asian","Volga-Ural","South Western Asian")
order_pop2 <- c("Northern Baltic","Southern Baltic","Eastern European Lowlands","Central European Lowlands",
                "Karelian","Eastern Finland","Western Finland", "Scandinavian", "Carpathian" ,
                "Dinaric Balkan","North Western Iberian","Italian")

metadata_mw$Group <- factor(metadata_mw$Group, levels = order_pop, ordered = T)
metadata_mw$GroupID <- factor(metadata_mw$GroupID, levels = order_popid, ordered = T)

groups$Group <- factor(groups$Group, levels = order_pop, ordered = T)
groups$GroupID <- factor(groups$GroupID, levels = order_popid, ordered = T)


metadata_mw2 <- merge(metadata_mw, groups[, c("Group","colour","SuperGroup")], by="Group")

metadata_mw2$SuperGroup <- factor(metadata_mw2$SuperGroup, levels=c("Fennoscandia","Eastern","Southern","Asian"), ordered = T)

# READ F4 RATIO FOR DOG ANCESTRY CORRELATIONS ----
load("Files/f4ratio_all.RData") 
res_f4$SampleID <- res_f4$pop3
metadata_mw_f4 <- merge(metadata_mw,res_f4[,c("SampleID","alpha","z")], by="SampleID")
metadata_mw2_f4 <- merge(metadata_mw2,res_f4[,c("SampleID","alpha","z")], by="SampleID")

metadata_mw2_f4_european <- metadata_mw2_f4[metadata_mw2_f4$SuperGroup!="Asian",]

summary(lm(metadata_mw2_f4_european$Heterozygosity~metadata_mw2_f4_european$alpha))
summary(lm(metadata_mw2_f4_european$Heterozygosity~metadata_mw2_f4_european$alpha+metadata_mw2_f4_european$Group))

g <- ggplot(metadata_mw2_f4_european, aes(alpha,Heterozygosity))
g <- g + geom_point(aes(color=colour)) + theme_classic()+geom_smooth(method="lm",formula=y~x)+
  scale_color_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  theme(legend.title = element_blank())

g

g1 <- ggplot(metadata_mw2_f4_european[which(metadata_mw2_f4_european$Group%in%
                                             c("Italian","Eastern European Lowlands","Southern Baltic")),], 
            aes(alpha, Heterozygosity, color=colour))
g1 <- g1 + geom_point(aes(color=colour)) + theme_classic()+geom_smooth(method="lm")+ 
 # facet_grid(.~Group)+
  scale_color_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  theme(legend.title = element_blank())


g1

dog <- g / g1 + plot_layout(heights = c(1, 1)) + plot_annotation(tag_levels = 'A')
dog

ggsave("Plots/Dog_het.pdf", height = 13, width = 13)
ggsave("Plots/FigS26.pdf", height = 13, width = 14)

library(broom)

results <- metadata_mw2_f4_european %>%
  group_by(Group) %>%
  group_modify(~ broom::tidy(lm(Heterozygosity ~ alpha, data = .x)))

results

summary(lm(metadata_mw2_f4_european[metadata_mw2_f4_european$Group=="Italian",]$Heterozygosity~
             metadata_mw2_f4_european[metadata_mw2_f4_european$Group=="Italian",]$alpha))

summary(lm(metadata_mw2_f4_european[metadata_mw2_f4_european$Group=="Eastern European Lowlands",]$Heterozygosity~
             metadata_mw2_f4_european[metadata_mw2_f4_european$Group=="Eastern European Lowlands",]$alpha))

summary(lm(metadata_mw2_f4_european[metadata_mw2_f4_european$Group=="Southern Baltic",]$Heterozygosity~
             metadata_mw2_f4_european[metadata_mw2_f4_european$Group=="Southern Baltic",]$alpha))


metadata_mw2_f4_european %>%
  group_by(Group) %>%
  summarise(
    mean_Heterozygosity = mean(Heterozygosity, na.rm = TRUE),
    sd_Heterozygosity   = sd(Heterozygosity, na.rm = TRUE),        # optional
    n_samples = n()                             # optional
  ) %>% arrange(desc(mean_Heterozygosity))




g <- ggplot(metadata_mw2_f4, aes(log(Heterozygosity)))
g <- g +geom_density()+
  scale_color_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) 

g


metadata_mw_f4$alpha0 <- metadata_mw_f4$alpha
metadata_mw_f4[metadata_mw_f4$alpha<0,]$alpha0 <- 0



metadata_mw_f4 %>%
  #filter(z > 3) %>%        
  group_by(Group) %>%
  summarise(
    mean_dog = mean(alpha0, na.rm = TRUE),
    sd_dog   = sd(alpha0, na.rm = TRUE),        # optional
    n_samples = n()                             # optional
  ) %>% arrange(desc(mean_dog))

#PLOT matrix of correlations between Froh Het and del allele counts----
library(ggplot2)
library(GGally)

df <- metadata_mw2_f4[, c("Heterozygosity", "FROH", 
                   "Total_High","Total_Moderate","Total_Low",
                   "Hom_High","Hom_Moderate","Hom_Low")]

# ---- custom upper panel: correlation with significance stars ----
upper_cor <- function(data, mapping, ...){
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  
  # get variable names
#  xlab <- as_label(mapping$x)
 # ylab <- as_label(mapping$y)
  
  ct <- suppressWarnings(cor.test(x, y))     #  Pearson by default
  r  <- unname(ct$estimate)
  p  <- ct$p.value
  stars <- symnum(p,
                  cutpoints = c(0, .001, .01, .05, 1),
                  symbols   = c("***", "**", "*", ""))
  lab <- sprintf("Corr:\n%.3f%s", r, stars)
  
  ggplot() +
    # y variable name (top)
   # annotate("text", x = 0.5, y = 0.9, label = ylab, fontface = "bold", size = 3.5) +
    # correlation in the middle
    annotate("text", x = 0.5, y = 0.5, label = lab, size = 4) +
    # x variable name (bottom)
  #  annotate("text", x = 0.5, y = 0.1, label = xlab, fontface = "bold", size = 3.5) +
    theme_void()
}

# ---- build the plot ----
p <- ggpairs(
  df,
  diag  = list(continuous = "densityDiag"),
  lower = list(continuous = wrap("points", alpha = 0.5, size = 1.4)),
  upper = list(continuous = upper_cor),
  progress = FALSE
) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text  = element_text(size = 8)
  )

p

ggsave("Plots/Correlations.pdf", height = 13, width = 13)
ggsave("Plots/FigS25.pdf", height = 13, width = 14)
