#Heterozygosity Wolves 
library(ggplot2)
library(readxl)
library(tidyr)
library(dplyr)
library("patchwork")
library(ggtext)
library(ggpubr)

setwd("~/Documents/OneDrive - University of Copenhagen/ResQ_MSCA/")

metadata_mw <- read_excel("Files/Wolves_metadata_Juny2026.xlsx", sheet = 1)
samples <- metadata_mw$SampleID

groups <- read_excel("Files/Wolves_metadata_Juny2026.xlsx", sheet = 2)

order_popid <- c("EASIA","SWASIA","ZOO","KAR",
                 "EFIN","WFIN", "SCAN", "NBAL",
                 "SBAL", "CEL","VOUR","CAR" ,
                 "DINBAL","NWIB","ITA")

order_pop <- c("Karelian","Eastern Finland",
               "Western Finland", "Scandinavian", "Northern Baltic",
               "Southern Baltic","Carpathian" ,"Central European Lowlands",
               "Dinaric Balkan","North Western Iberian","Italian",
               "Eastern Asian","Volga-Ural","South Western Asian","Zoo")
order_pop2 <- c("Northern Baltic","Southern Baltic","Central European Lowlands","Volga-Ural",
                "Karelian","Eastern Finland","Western Finland", "Scandinavian","Zoo", "Carpathian" ,
                "Dinaric Balkan","North Western Iberian","Italian")

metadata_mw$Group <- factor(metadata_mw$Group, levels = order_pop, ordered = T)
metadata_mw$GroupID <- factor(metadata_mw$GroupID, levels = order_popid, ordered = T)

groups$Group <- factor(groups$Group, levels = order_pop, ordered = T)
groups$GroupID <- factor(groups$GroupID, levels = order_popid, ordered = T)


metadata_mw2 <- merge(metadata_mw, groups[, c("Group","colour","SuperGroup")], by="Group")

metadata_mw2$SuperGroup <- factor(metadata_mw2$SuperGroup, levels=c("Eastern","Fennoscandia","Southern","Asian","Zoo"), ordered = T)


metadata_mw2_rm <- metadata_mw2 %>%
  filter(is.na(Remove) | Remove!="Yes")


f4_ratio<- read.delim("Files/f4ratiostats_all.txt", sep="\t", header = F)

hybrids <- f4_ratio[f4_ratio$V14=="Yes",]$V3

# Some stats on heterozygosity----
summary <- metadata_mw2 %>%
  filter(is.na(Remove) | Remove!="Yes")%>%
  group_by(Group) %>%
  summarise(
    mean_Heterozygosity = mean(Het_VCF, na.rm = TRUE),
    sd_Heterozygosity   = sd(Het_VCF, na.rm = TRUE),        
    mean_Heterozygosity_angsd = mean(Het_BAM, na.rm = TRUE),
    mean_Froh   = mean(FROH, na.rm = TRUE),     
    mean_ExtraROH = mean(ExtraLong, na.rm = TRUE), 
    mean_LongROH = mean(Long, na.rm = TRUE), 
    mean_InterROH = mean(Intermediate, na.rm = TRUE), 
    mean_ShortROH = mean(Short, na.rm = TRUE), 
    mean_COVERAGE = mean(Coverage, na.rm = TRUE),
    n_samples = n()                        
  ) %>% arrange(desc(mean_Heterozygosity))

summary_supergroups <- metadata_mw2 %>%
  filter(is.na(Remove) | Remove!="Yes")%>%
  group_by(SuperGroup) %>%
  summarise(
    mean_Heterozygosity = mean(Het_VCF, na.rm = TRUE),
    sd_Heterozygosity   = sd(Het_VCF, na.rm = TRUE),        
    mean_Heterozygosity_angsd = mean(Het_BAM, na.rm = TRUE),
    mean_Froh   = mean(FROH, na.rm = TRUE),     
    mean_ExtraROH = mean(ExtraLong, na.rm = TRUE), 
    mean_LongROH = mean(Long, na.rm = TRUE), 
    mean_InterROH = mean(Intermediate, na.rm = TRUE), 
    mean_ShortROH = mean(Short, na.rm = TRUE), 
    mean_COVERAGE = mean(Coverage, na.rm = TRUE),
    n_samples = n()                        
  ) %>% arrange(desc(mean_Heterozygosity))

View(summary)
average_values <- metadata_mw2 %>%
  filter(is.na(Remove) | Remove!="Yes")%>%
  summarise(
    mean_Heterozygosity = mean(Het_VCF, na.rm = TRUE),
    mean_Heterozygosity_angsd = mean(Het_BAM, na.rm = TRUE),
    mean_Froh   = mean(FROH, na.rm = TRUE),     
    n_samples = n()                        
  ) %>% arrange(desc(mean_Heterozygosity))
average_values


# Plot heterozygosity per group without ASIAN ----
g<- ggplot(metadata_mw2_rm[metadata_mw2_rm$SuperGroup!="Asian",], aes(Group, Het_VCF, color=colour))
g <- g + geom_point(position = "jitter", alpha=0.6)+
  geom_boxplot(outlier.shape = NA, fill=NA )+
  facet_grid(.~SuperGroup, scales="free", space="free")+
  theme_classic()+
  scale_color_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  xlab("")+ ylab("Relative Heterozygosity")+
  theme(axis.text.x = element_text(angle=45, hjust=1), 
        legend.position = "none")
g

ggsave("Plots/Heterozygosity_Groups.pdf", width = 15, height = 7)


# plot heterozygosity per group including ASIAN -----
g<- ggplot(metadata_mw2_rm, aes(Group, Het_VCF, color=colour,group=colour, shape=Status))
g <- g + geom_point(position = "jitter", alpha=0.6)+
  geom_boxplot(outlier.shape = NA, fill=NA )+
  geom_hline(yintercept = average_values$mean_Heterozygosity, linetype="dashed", color="grey")+
  theme_classic()+ 
  facet_grid(.~SuperGroup, scales="free", space="free")+
  scale_color_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  scale_shape_manual(values=c(16,17))+
  xlab("")+ ylab("Relative Heterozygosity")+
  guides(color=guide_legend(title="Population"))+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none")
g

ggsave("Plots/Heterozygosity_all.pdf", width = 10, height = 4)


# plot heterozygosity ANGSD per group including ASIAN -----
g<- ggplot(metadata_mw2_rm, aes(Group, Het_BAM, color=colour,group=colour, shape=Status, label=SampleID))
g <- g + geom_point(position = "jitter", alpha=0.6)+
  geom_text(data=metadata_mw2[metadata_mw2$Het_BAM>0.004,], alpha=0.6)+
  geom_boxplot(outlier.shape = NA, fill=NA )+
  geom_hline(yintercept = average_values$mean_Heterozygosity_angsd, linetype="dashed", color="grey")+
  theme_classic()+ 
  facet_grid(.~SuperGroup, scales="free", space="free")+
  scale_color_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  scale_shape_manual(values=c(16,17))+
  xlab("")+ ylab("Heterozygosity ANGSD")+
  guides(color=guide_legend(title="Population"))+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none")
g

ggsave("Plots/Heterozygosity_ANGSD_all.pdf", width = 10, height = 4)

col_vec_as <- setNames(groups$colour, groups$GroupID)
xlab_vec_as <- setNames(
  paste0("<span style='color:", col_vec_as, "'>", groups$GroupID, "</span>"),
  groups$GroupID
)

g_angsd<- ggplot(metadata_mw2_rm[metadata_mw2_rm$Het_BAM<0.004,], aes(GroupID, Het_BAM, color=colour,group=colour, shape=Status))
g_angsd <- g_angsd + geom_point(position = "jitter", alpha=0.6)+
  geom_boxplot(outlier.shape = NA, fill=NA )+
  geom_hline(yintercept = average_values$mean_Heterozygosity_angsd, linetype="dashed", color="grey")+
  theme_classic()+ 
  facet_grid(.~SuperGroup, scales="free", space="free")+
  scale_color_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  scale_shape_manual(values=c(16,17))+
  xlab("")+ ylab("Heterozygosity ANGSD")+
  scale_x_discrete(labels = xlab_vec_as)+
  guides(color=guide_legend(title="Population"))+
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none")
g_angsd

ggsave("Plots/Heterozygosity_ANGSD_nooutliers.pdf", width = 10, height = 4)


# with weird outliers

fit <- lm(Het_VCF ~ Het_BAM, data = metadata_mw2)
summary(fit)

g_cor<- ggplot(metadata_mw2, aes(Het_VCF, Het_BAM, label=SampleID))
g_cor <- g_cor + geom_point(position = "jitter", alpha=0.6, aes(color=colour,group=colour, shape=Status))+
  geom_text(data=metadata_mw2[metadata_mw2$Het_BAM>0.002,], alpha=0.6)+
  geom_text(data=metadata_mw2[metadata_mw2$Het_VCF>0.11,], alpha=0.6)+
  geom_smooth(method = "lm", se = TRUE,color="black") +
  ylab("Heterozygosity ANGSD") + xlab("Relative Heterozygosity VCF")+
  theme_classic()+ 
  #facet_grid(.~SuperGroup, scales="free", space="free")+
  scale_color_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  scale_shape_manual(values=c(16,17))+
  guides(color=guide_legend(title="Population"))
g_cor

ggsave("Plots/Heterozygosity_ANGSD_VCF.pdf", width = 10, height = 7)

fit <- lm(Coverage ~ Het_BAM, data = metadata_mw2[metadata_mw2$Het_BAM<0.002,])
summary(fit)


g_cov<- ggplot(metadata_mw2, aes(Coverage, Het_BAM))
g_cov <- g_cov + geom_point(aes(color=colour,group=colour, shape=Status, label=SampleID),position = "jitter", alpha=0.6)+
  geom_text(data=metadata_mw2[metadata_mw2$Het_BAM>0.002,], aes(color=colour,group=colour, shape=Status, label=SampleID), alpha=0.6)+
  theme_classic()+ 
  #facet_grid(.~SuperGroup, scales="free", space="free")+
  scale_color_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  scale_shape_manual(values=c(16,17))+ ylab("Heterozygosity ANGSD")+
  geom_smooth(data=metadata_mw2[metadata_mw2$Het_BAM<0.002,], method = "lm", se = TRUE,color="black") +
  theme(legend.position = "none")+
  guides(color=guide_legend(title="Population"))
g_cov

g_angsd/g_cov

het_angsd_all <- g_angsd / g_cov + plot_layout(heights = c(1, 2)) + plot_annotation(tag_levels = 'A')
ggsave("Plots/Het_ANGSD_Cov_SFig.pdf", width = 7, height = 8)

g<- ggplot(metadata_mw2, aes(Coverage, Het_BAM, color=Coverage, label=SampleID))
g <- g + geom_point(position = "jitter", alpha=0.6)+
  geom_text(data=metadata_mw2[metadata_mw2$Het_BAM>0.002,], alpha=0.6)+
  theme_classic()+ 
  #facet_grid(.~SuperGroup, scales="free", space="free")+
  scale_color_viridis_c()
g


g<- ggplot(metadata_mw2[metadata_mw2$Coverage<4,], aes(Coverage, Het_BAM, color=Coverage, label=SampleID))
g <- g + geom_point(position = "jitter", alpha=0.6)+
  geom_text(alpha=0.6)+
  theme_classic()+ 
  scale_color_viridis_c()
g


#no transitions
g_angsd<- ggplot(metadata_mw2_rm[which(metadata_mw2_rm$Het_BAM_notrans<0.004),], aes(GroupID, Het_BAM_notrans, color=colour,group=colour, shape=Status))
g_angsd <- g_angsd + geom_point(position = "jitter", alpha=0.6)+
  geom_boxplot(outlier.shape = NA, fill=NA )+
  geom_hline(yintercept = average_values$mean_Heterozygosity_angsd, linetype="dashed", color="grey")+
  theme_classic()+ 
  facet_grid(.~SuperGroup, scales="free", space="free")+
  scale_color_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  scale_shape_manual(values=c(16,17))+
  xlab("")+ ylab("Heterozygosity ANGSD (no transitions)")+
  scale_x_discrete(labels = xlab_vec_as)+
  guides(color=guide_legend(title="Population"))+
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none")
g_angsd

ggsave("Plots/Heterozygosity_ANGSD_nooutliers_notrans.pdf", width = 10, height = 4)

g_cov<- ggplot(metadata_mw2, aes(Coverage, Het_BAM_notrans))
g_cov <- g_cov + geom_point(aes(color=colour,group=colour, shape=Status, label=SampleID),position = "jitter", alpha=0.6)+
  geom_text(data=metadata_mw2[metadata_mw2$Het_BAM>0.002,], aes(color=colour,group=colour, shape=Status, label=SampleID), alpha=0.6)+
  theme_classic()+ 
  #facet_grid(.~SuperGroup, scales="free", space="free")+
  scale_color_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  scale_shape_manual(values=c(16,17))+ ylab("Heterozygosity ANGSD (no transitions)")+
  geom_smooth(data=metadata_mw2[metadata_mw2$Het_BAM<0.002,], method = "lm", se = TRUE,color="black") +
  theme(legend.position = "none")+
  guides(color=guide_legend(title="Population"))
g_cov

g_angsd/g_cov

het_angsd_all <- g_angsd / g_cov + plot_layout(heights = c(1, 2)) + plot_annotation(tag_levels = 'A')
ggsave("Plots/Het_ANGSD_Cov_notrans.pdf", width = 7, height = 8)



# plot heterozygosity per group including ASIAN -----
col_vec_as <- setNames(groups$colour, groups$GroupID)
xlab_vec_as <- setNames(
  paste0("<span style='color:", col_vec_as, "'>", groups$GroupID, "</span>"),
  groups$GroupID
)

g<- ggplot(metadata_mw2_rm, aes(GroupID, Het_VCF, color=colour,group=colour, shape=Status))
g <- g + geom_point(position = "jitter", alpha=0.6)+
  geom_boxplot(outlier.shape = NA, fill=NA )+
  geom_hline(yintercept = average_values$mean_Heterozygosity, linetype="dashed", color="grey")+
  theme_classic()+ 
  facet_grid(.~SuperGroup, scales="free", space="free")+
  scale_color_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  scale_shape_manual(values=c(16,17))+
  scale_x_discrete(labels = xlab_vec_as)+
  xlab("")+ ylab("Relative Heterozygosity")+
  guides(color=guide_legend(title="Population"))+
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none")
g


# plot FROH per group
g1<- ggplot(metadata_mw2_rm, aes(GroupID, FROH, color=colour, group=colour,shape=Status))
g1 <- g1 + geom_point(position = "jitter", alpha=0.6)+
  geom_boxplot(outlier.shape = NA, fill=NA )+
  geom_hline(yintercept = average_values$mean_Froh, linetype="dashed", color="grey")+
  theme_classic()+ 
  facet_grid(.~SuperGroup, scales="free", space="free")+
  scale_color_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  scale_shape_manual(values=c(16,17))+
  scale_x_discrete(labels = xlab_vec_as)+
  xlab("")+ ylab(expression(F[ROH]))+
  guides(color=guide_legend(title="Population"))+
theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
                 legend.position = "none")
g1

ggsave("Plots/FROH_all.pdf", width = 10, height = 4)


het_froh_all <- g / g1 + plot_layout(heights = c(1, 1)) + plot_annotation(tag_levels = 'A')
het_froh_all
ggsave("Plots/FigureS20.pdf", width = 10, height = 8)


# Heterozygosity excluding hybrids
f4_ratio$SampleID <- f4_ratio$V3
metadata_mw2_rm_hyb <- merge(metadata_mw2_rm,f4_ratio[,c("SampleID","V14")], by="SampleID")




g<- ggplot(metadata_mw2_rm[-which(metadata_mw2_rm$SampleID%in%hybrids),], aes(GroupID, Het_VCF, color=colour,group=colour, shape=Status))
g <- g + geom_point(position = "jitter", alpha=0.6)+
  geom_boxplot(outlier.shape = NA, fill=NA )+
  geom_hline(yintercept = average_values$mean_Heterozygosity, linetype="dashed", color="grey")+
  theme_classic()+ 
  facet_grid(.~SuperGroup, scales="free", space="free")+
  scale_color_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  scale_shape_manual(values=c(16,17))+
  scale_x_discrete(labels = xlab_vec_as)+
  xlab("")+ ylab("Relative Heterozygosity")+
  guides(color=guide_legend(title="Population"))+
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none")
g





g<- ggplot(metadata_mw2_rm_hyb, aes(GroupID, Het_VCF, color=colour,group=colour, shape=V14))
g <- g +   geom_point(aes(group = interaction(GroupID, V14)),
                      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
                      alpha = 0.6) +
 # geom_boxplot(outlier.shape = NA, fill=NA )+
  geom_boxplot(aes(group = interaction(GroupID, V14)), outlier.shape = NA, fill = NA,
               position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = average_values$mean_Heterozygosity, linetype="dashed", color="grey")+
  stat_compare_means(aes(group = V14), method = "wilcox.test", label = "p.signif") +
  theme_classic()+ 
  facet_grid(.~SuperGroup, scales="free", space="free")+
  scale_color_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  scale_shape_manual(values=c(16,17))+
  xlab("")+ ylab("Relative Heterozygosity")+
  guides(color=guide_legend(title="Population"))+
  scale_x_discrete(labels = xlab_vec_as)+
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none")
g
ggsave("Plots/FigureS35.pdf", width = 10, height = 4)

