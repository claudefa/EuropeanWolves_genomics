#Heterozygosity Wolves RESQ
library(ggplot2)
library(readxl)
library(tidyr)
library(dplyr)
library("patchwork")

setwd("~/Documents/OneDrive - University of Copenhagen/ResQ_MSCA/")

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


# Some stats on heterozygosity----
metadata_mw2 %>%
  group_by(Group) %>%
  summarise(
    mean_Heterozygosity = mean(Heterozygosity, na.rm = TRUE),
    sd_Heterozygosity   = sd(Heterozygosity, na.rm = TRUE),        # optional
    n_samples = n()                             # optional
  ) %>% arrange(desc(mean_Heterozygosity))

average_values <- metadata_mw2 %>%
  summarise(
    mean_Heterozygosity = mean(Heterozygosity, na.rm = TRUE),
    mean_Froh   = sd(FROH, na.rm = TRUE),        # optional
    n_samples = n()                             # optional
  ) %>% arrange(desc(mean_Heterozygosity))
average_values


# Plot heterozygosity per group without ASIAN ----
g<- ggplot(metadata_mw2[metadata_mw2$SuperGroup!="Asian",], aes(Group, Heterozygosity, color=colour))
g <- g + geom_point(position = "jitter", alpha=0.6)+
  geom_boxplot(outlier.shape = NA, fill=NA )+
  facet_grid(.~SuperGroup, scales="free", space="free")+
  theme_classic()+
  scale_color_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  xlab("")+ ylab("Relative Heterozygosity")+
  theme(axis.text.x = element_text(angle=45, hjust=1), 
        legend.position = "none")
g

ggsave("Plots/Heterozygosity_Groups_52M.pdf", width = 15, height = 7)


# plot heterozygosity per group including ASIAN -----
g<- ggplot(metadata_mw2, aes(Group, Heterozygosity, color=colour,group=colour, shape=Status))
g <- g + geom_point(position = "jitter", alpha=0.6)+
  geom_boxplot(outlier.shape = NA, fill=NA )+
  geom_hline(yintercept = average_values$mean_Heterozygosity, linetype="dashed", color="grey")+
  theme_classic()+ 
  facet_grid(.~SuperGroup, scales="free", space="free")+
  scale_color_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  scale_shape_manual(values=c(17,16))+
  xlab("")+ ylab("Relative Heterozygosity")+
  guides(color=guide_legend(title="Population"))+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none")
g

ggsave("Plots/Heterozygosity_all_52M.pdf", width = 10, height = 4)


# plot FROH per group
g1<- ggplot(metadata_mw2, aes(Group, FROH, color=colour, group=colour,shape=Status))
g1 <- g1 + geom_point(position = "jitter", alpha=0.6)+
  geom_boxplot(outlier.shape = NA, fill=NA )+
  geom_hline(yintercept = average_values$mean_Froh, linetype="dashed", color="grey")+
  theme_classic()+ 
  facet_grid(.~SuperGroup, scales="free", space="free")+
  scale_color_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  scale_shape_manual(values=c(17,16))+
  xlab("")+ ylab(expression(F[ROH]))+
  guides(color=guide_legend(title="Population"))+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none")
g1

ggsave("Plots/Heterozygosity_all_52M.pdf", width = 10, height = 4)


het_froh_all <- g / g1 + plot_layout(heights = c(1, 1)) + plot_annotation(tag_levels = 'A')
het_froh_all
ggsave("Plots/Heterozygosity_Froh_52M_all.pdf", width = 10, height = 8)
ggsave("Plots/FigureS20.pdf", width = 10, height = 8)







