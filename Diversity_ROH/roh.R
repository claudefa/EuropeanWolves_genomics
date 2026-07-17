# ROH modern wolves
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
library(cowplot)

setwd("~/Documents/OneDrive - University of Copenhagen/ResQ_MSCA/")

# Metadata ----
metadata_mw <- read_excel("Files/Wolves_metadata_Juny2026.xlsx", sheet = 1)
metadata_mw<-metadata_mw %>%
  filter(is.na(Remove) | Remove!="Yes")

samples <- metadata_mw$SampleID

groups <- read_excel("Files/Wolves_metadata_Juny2026.xlsx", sheet = 2)

order_popid <- c("EASIA","VOUR","SWASIA","ZOO","KAR",
                 "EFIN","WFIN", "SCAN", "NBAL",
                 "SBAL","CAR" ,"CEL",
                 "DINBAL","NWIB","ITA")

order_pop <- c("Karelian","Eastern Finland",
               "Western Finland", "Scandinavian", "Northern Baltic",
               "Southern Baltic","Carpathian" ,"Central European Lowlands","Volga-Ural",
               "Dinaric Balkan","North Western Iberian","Italian",
               "Eastern Asian","South Western Asian","Zoo")
order_pop2 <- c("Northern Baltic","Southern Baltic","Central European Lowlands",
                "Karelian","Eastern Finland","Western Finland", "Scandinavian","Zoo", "Carpathian" ,
                "Dinaric Balkan","North Western Iberian","Italian")

metadata_mw$Group <- factor(metadata_mw$Group, levels = order_pop, ordered = T)
metadata_mw$GroupID <- factor(metadata_mw$GroupID, levels = order_popid, ordered = T)

groups$Group <- factor(groups$Group, levels = order_pop, ordered = T)
groups$GroupID <- factor(groups$GroupID, levels = order_popid, ordered = T)


metadata_mw2 <- merge(metadata_mw, groups[, c("Group","colour","SuperGroup")], by="Group")
metadata_mw2$SuperGroup <- factor(metadata_mw2$SuperGroup, levels=c("Fennoscandia","Eastern","Southern","Asian","Zoo"), ordered = T)

# ROH ------
roh_file <- read.table("Files/ROH/modernWolf_AF_map.txt") # with estimated AF from the VCF and Dog recomb map

colnames(roh_file) <- c("V1","SampleID", "CHR","START","END","length","markers", "quality")
assembly <- read_excel("Files/Wolves_metadata_Juny2026.xlsx", sheet = 3)
chr_order <- assembly$CHR
assembly$CHR <- factor( assembly$CHR, ordered = T, levels = chr_order)

# Quality and size evaluation -----

g <- ggplot(roh_file[roh_file$SampleID=="W10",], aes(quality,length))
g <- g+ geom_point() + theme_light()+
  geom_hline(yintercept = 500000, linetype="dashed", color="red") +
  geom_vline(xintercept = 50, linetype="dashed",color="red") +
  xlim(c(50,100))
g

#minimum filtering to remove the untrustworthy rohs ----
roh_final <- roh_file[roh_file$quality > 50&roh_file$length>=500000,]# keep only rohs > 500000 and quality >50
roh_final <- roh_final[which(roh_final$CHR%in%chr_order),] # only autosomes

# Add age of ROH -----
GenTime=4.5 #Smeds et al 2019
r=1 # 0.97 dog genome   (Campbell et al., 2016) (Pacheco et al. 2022)
# formula is L = 100/2t cM (from Thompson 2013, Genetics 194, 301-326)
roh_final$Tgen=(100/((roh_final$length/1000000)*r))/2
roh_final$Tyears=roh_final$Tgen*GenTime


## Percentage of the genome in ROH per sample ----
sum(assembly$SIZE) # size of the genome

roh_perc <- lapply(1:length(samples), function(i) {
  df2 <- roh_final[which(roh_final$SampleID == samples[i]),]
  data.frame(SampleID=samples[i],
             Total=sum(df2$length/2203764842),
             Short=sum(df2[df2$length>=500000 & df2$length<=1000000,]$length/2203764842), 
             Intermediate=sum(df2[df2$length>1000000 & df2$length<=2500000,]$length/2203764842), 
             Long=sum(df2[df2$length>2500000 & df2$length<=50000000,]$length/2203764842), 
             ExtraLong=sum(df2[df2$length>5000000,]$length/2203764842))
})


roh_perc_df <- do.call(rbind, roh_perc)

#add metadata info
roh_perc_df_m <- merge(roh_perc_df, metadata_mw2[,c("SampleID","Location","Group","GroupID","SuperGroup")], by="SampleID")
roh_perc_df_m$Group <- factor(roh_perc_df_m$Group, levels = order_pop, ordered = T)
roh_perc_df_m$GroupID <- factor(roh_perc_df_m$GroupID, levels = order_popid, ordered = T)

#Statistics on size of ROH ---
mean(roh_perc_df_m[roh_perc_df_m$SuperGroup=="Fennoscandia",]$Short)*100
mean(roh_perc_df_m[roh_perc_df_m$SuperGroup=="Fennoscandia",]$ExtraLong)*100 

mean(roh_perc_df_m[roh_perc_df_m$Group=="Scandinavian",]$Short)*100
mean(roh_perc_df_m[roh_perc_df_m$Group=="Scandinavian",]$ExtraLong)*100 

mean(roh_perc_df_m[roh_perc_df_m$Group=="Italian",]$Short)*100
mean(roh_perc_df_m[roh_perc_df_m$Group=="Italian",]$ExtraLong)*100 


mean(roh_perc_df_m[roh_perc_df_m$Group=="Central European Lowlands",]$Short)*100
mean(roh_perc_df_m[roh_perc_df_m$Group=="Central European Lowlands",]$ExtraLong)*100 

# Save FROH----
write.csv2(roh_perc_df_m, "Files/Froh.csv", quote = F, row.names = F)

#Formatting for plotting ----
roh_perc_df_gather <- gather(roh_perc_df_m, key, value, -Total, -SampleID, -Location, -Group,-GroupID, -SuperGroup)

roh_perc_df_gather$key <- gsub("Short","500Kb-1Mb", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("Intermediate","1-2.5Mb", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("ExtraLong",">5Mb", roh_perc_df_gather$key )
roh_perc_df_gather$key <- gsub("Long","2.5-5Mb", roh_perc_df_gather$key )

order <- c("500Kb-1Mb","1-2.5Mb", "2.5-5Mb",">5Mb")

roh_perc_df_gather$key <- factor(roh_perc_df_gather$key, levels = order, ordered = TRUE)


# Individual ROH plot
pa <- ggplot(roh_perc_df_gather, aes(SampleID,value*100,fill=key))
pa <- pa + geom_col()+ theme_bw() +  xlab("")+ 
  ylab("% of Genome in RoH") + facet_grid(.~Group, scales="free", space="free",switch="both")+
  theme(legend.position = "right", strip.background = element_blank(), axis.text.y = element_text(hjust=0.5)) + 
  scale_fill_manual(values=rev(c("#fcd021","#193563","#8cb4f5", "#cce0ff")), name="RoH Size") +
  scale_y_continuous(limits = c(0,60), expand = c(0, 0)) +
  theme(axis.text.x = element_blank(), axis.line.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank(), strip.text = element_text(angle=90, hjust = 1))
pa
ggsave("Plots/Summary_rohs_modernwolves_AF_recombmap.pdf", width = 25, height = 5)


# Plot by group FROH
roh_perc_df_m2 <- merge(roh_perc_df_m, groups[, c("Group","colour")], by="Group")

froh<- ggplot(roh_perc_df_m2[roh_perc_df_m2$SuperGroup!="Asian",], aes(Group, Total*100, color=colour))
froh <- froh + geom_point(position = "jitter", alpha=0.6)+
  geom_boxplot(outlier.shape = NA, fill=NA )+
  facet_grid(.~SuperGroup, scales="free", space="free")+
  theme_classic()+ scale_shape_manual(values=c(2,1,3))+
  scale_color_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  xlab("")+ ylab("% of Genome in ROH")+
  guides(color=guide_legend(title="Population"))+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none")
froh

# ROH by size ------
roh_perc_df_gather_m <- merge(roh_perc_df_gather, groups[, c("Group","colour")], by="Group")

g<- ggplot(roh_perc_df_gather_m, aes(Group, value*100, color=colour, group=interaction(Group)))
g <- g + geom_point(position = "jitter", alpha=0.6)+
  geom_boxplot(outlier.shape = NA, fill=NA )+  facet_wrap(.~key, ncol=4)+
  theme_classic()+ scale_shape_manual(values=c(2,1,3))+
  scale_color_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  xlab("")+ ylab("% of Genome in ROH")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none")
g

ggsave("Plots/Total_rohs_bysize_modernwolves_AF_recombmap.pdf", width = 16, height = 7)


roh_perc_df_gather_m$Group <- factor(roh_perc_df_gather_m$Group, levels = order_pop, ordered = TRUE)

g_size<- ggplot(roh_perc_df_gather_m[roh_perc_df_gather_m$SuperGroup!="Asian",],
           aes(Group, value*100, color=colour))
g_size <- g_size + geom_point(position = "jitter", alpha=0.6)+
  geom_boxplot(outlier.shape = NA, fill=NA )+  facet_wrap(.~key, ncol=4, scales="free_y")+
  theme_classic()+
  scale_color_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  xlab("")+ ylab("% of Genome in ROH")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none")
g_size
ggsave("Plots/Total_rohs_bysize_modernwolves_AF_recombmap_noAsian.pdf", width = 16, height = 7)


# ROH by time ------
roh_time <- lapply(1:length(samples), function(i) {
  df2 <- roh_final[which(roh_final$SampleID == samples[i]),]
  data.frame(SampleID=samples[i],
             cat2=sum(df2[df2$Tgen>=0 & df2$Tgen<10,]$length/2203764842), 
             cat3=sum(df2[df2$Tgen>=10 & df2$Tgen<20,]$length/2203764842), 
             cat4=sum(df2[df2$Tgen>=20&df2$Tgen<50,]$length/2203764842),
             cat5=sum(df2[df2$Tgen>=50,]$length/2203764842))
})


roh_time_df <- do.call(rbind, roh_time)

#add metadata info
roh_time_df_m <- merge(roh_time_df, metadata_mw2[,c("SampleID","Location","Group","GroupID","SuperGroup", "Status")], by="SampleID")
roh_time_df_m$Group <- factor(roh_time_df_m$Group, levels = order_pop, ordered = T)
roh_time_df_m$GroupID <- factor(roh_time_df_m$GroupID, levels = order_popid, ordered = T)

roh_time_df_gather <- gather(roh_time_df_m, key, value, -SampleID, -Location, -Group,-GroupID, -SuperGroup, -Status)

roh_time_df_gather$key <- gsub("cat2","<10 Gen (> 5Mb)", roh_time_df_gather$key )
roh_time_df_gather$key <- gsub("cat3",">10 - 20 Gen (2.5-5Mb)", roh_time_df_gather$key )
roh_time_df_gather$key <- gsub("cat4",">20 - 50 Gen (1-2.5Mb)", roh_time_df_gather$key )
roh_time_df_gather$key <- gsub("cat5",">50 Gen (500-1Mb)", roh_time_df_gather$key )

order <- c("<10 Gen (> 5Mb)", ">10 - 20 Gen (2.5-5Mb)",">20 - 50 Gen (1-2.5Mb)",">50 Gen (500-1Mb)")

roh_time_df_gather$key <- factor(roh_time_df_gather$key, levels = order, ordered = TRUE)
roh_time_df_gather_m <- merge(roh_time_df_gather, groups[, c("Group","colour")], by="Group")

roh_time_df_gather_m$Group <- factor(roh_time_df_gather_m$Group, levels = order_pop, ordered = TRUE)


g_age<- ggplot(roh_time_df_gather_m[roh_time_df_gather_m$SuperGroup!="Asian",],
           aes(Group, value*100, color=colour))
g_age <- g_age + geom_point(position = "jitter", alpha=0.6)+
  geom_boxplot(outlier.shape = NA, fill=NA )+  facet_wrap(.~key,scales="free_y", nrow=1)+
  theme_classic()+ 
  scale_color_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  xlab("")+ ylab("% of Genome in ROH")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none")
g_age

ggsave("Plots/Total_rohs_byage_modernwolves_AF_recombmap_noAsian.pdf", width = 16, height = 7)

