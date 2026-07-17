# Genetic Load with SNPEff ---
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(readxl)
library(bootstrap)
library(ggtext)
library(patchwork)
library(data.table)

setwd("~/Documents/OneDrive - University of Copenhagen/ResQ_MSCA/")

#Metadata
metadata_mw <- read_excel("Files/Wolves_metadata_Juny2026.xlsx", sheet = 1)
metadata_mw$Remove<- as.character(metadata_mw$Remove)
metadata_mw <- metadata_mw[-which(metadata_mw$Remove=="Yes"),]
samples <- metadata_mw$SampleID

groups <- read_excel("Files/Wolves_metadata_Juny2026.xlsx", sheet = 2)

order_popid <- c("EASIA","SWASIA","ZOO","KAR",
                 "EFIN","WFIN", "SCAN", "NBAL",
                 "SBAL","CEL","VOUR","CAR" ,
                 "DINBAL","NWIB","ITA")

order_pop <- c("Karelian","Eastern Finland",
               "Western Finland", "Scandinavian", "Northern Baltic",
               "Southern Baltic","Carpathian" ,"Central European Lowlands","Volga-Ural",
               "Dinaric Balkan","North Western Iberian","Italian",
               "Eastern Asian","South Western Asian","Zoo")
order_pop2 <- c("Northern Baltic","Southern Baltic","Central European Lowlands","Volga-Ural",
                "Karelian","Eastern Finland","Western Finland", "Scandinavian","Zoo", "Carpathian" ,
                "Dinaric Balkan","North Western Iberian","Italian")

metadata_mw$Group <- factor(metadata_mw$Group, levels = order_pop, ordered = T)
metadata_mw$GroupID <- factor(metadata_mw$GroupID, levels = order_popid, ordered = T)

groups$Group <- factor(groups$Group, levels = order_pop, ordered = T)
groups$GroupID <- factor(groups$GroupID, levels = order_popid, ordered = T)


metadata_mw2 <- merge(metadata_mw, groups[, c("Group","colour","SuperGroup")], by="Group")

metadata_mw2$SuperGroup <- factor(metadata_mw2$SuperGroup, levels=c("Eastern","Fennoscandia","Southern","Asian","Zoo"), ordered = T)


outgroups_id <- c("AlabamaCoyote","Alaska2","Algerian_grey_wolf","AndeanFox","BerlinZoo",
               "C_adustus","C_mesomelas","CaliforniaCoyote","Ethiopian_grey_wolf",
               "KenyaJackal","Kruger","Mexico","MW004","MW613","MW615","MW616",
               "MW617","MW618","MW621","MW670","MW672","MW673","MW674","MW675",
               "MW676","MW677","MW678","MW689","MW690","MW691","MW692","Syrian_jackal",
               "Zimbabwe1")

outgroup_g <- c("Coyote","Coyote","African-golden-wolf","Andean-fox","Dhole",
                "Side-Striped-Jackal","Black-Backed-Jackal","Coyote",
                "African-golden-wolf","African-golden-wolf","African-hunting-dog",
                "Coyote","Golden-jackal","Golden-jackal","Golden-jackal",
                "Golden-jackal","Golden-jackal","Golden-jackal","Golden-jackal",
                "Golden-jackal","Golden-jackal","Golden-jackal","Golden-jackal",
                "Golden-jackal","Golden-jackal","Golden-jackal","Golden-jackal",
                "Golden-jackal","Golden-jackal","Golden-jackal","Golden-jackal",
                "Golden-jackal","African-hunting-dog")


outgroup <- cbind.data.frame(outgroups_id,outgroup_g)  

samples_onlyWolves <- metadata_mw$SampleID

metadata_mw2$LON <- as.numeric(metadata_mw2$Longitude)
metadata_mw2$LAT <- as.numeric(metadata_mw2$Latitude)

groups_noAsian <- groups[groups$Group %in% order_pop2, ]
groups_noAsian <- groups_noAsian[match(order_pop2, groups_noAsian$Group), ]

groups <- groups[match(order_pop, groups$Group), ]

# Obtain counts for all groups ------
# read files for each category and remove sites with all sites have 0|0
samples<- read.table("Files/GenLoad/SNPeff/Samples_1117")

high <- read.table("Files/GenLoad/SNPeff/high.txt")
colnames(high) <- c("Chrom","Pos","REF","ALT", samples$V1)
high$chrompos<- paste0(high$Chrom,":",high$Pos)
keep <- apply(high[5:1121], 1, function(x) length(unique(x[!is.na(x)])) != 1)
high_nonvariable <- high[keep, ]
high_nonvariable$Type <- "High"

moderate <- read.table("Files/GenLoad/SNPeff/moderate.txt")
colnames(moderate) <- c("Chrom","Pos","REF","ALT", samples$V1)
moderate$chrompos<- paste0(moderate$Chrom,":",moderate$Pos)
keep <- apply(moderate[5:1121], 1, function(x) length(unique(x[!is.na(x)])) != 1)
moderate_nonvariable <- moderate[keep, ]
moderate_nonvariable$Type <- "Moderate"

low <- read.table("Files/GenLoad/SNPeff/low.txt")
colnames(low) <- c("Chrom","Pos","REF","ALT", samples$V1)
low$chrompos<- paste0(low$Chrom,":",low$Pos)
keep <- apply(low[5:1121], 1, function(x) length(unique(x[!is.na(x)])) != 1)
low_nonvariable <- low[keep, ]
low_nonvariable$Type <- "Low"

# Put all categories in a single dataframe
total_file <- rbind.data.frame(high_nonvariable,moderate_nonvariable,low_nonvariable)
total_sites_wononvariable <- length(total_file$Chrom)

# Function to add ancestral state ----
add_ref_ancestral_filter <- function(df, outgroups, min_support = 2) {
  # Convert genotypes to allele 0/1 or NA
  get_allele <- function(g) {
    ifelse(g == "0|0", "0",
           ifelse(g == "1|1", "1", NA))
  }
  
  # Extract 0/1 allele calls for tier1 outgroups
  og <- as.data.frame(lapply(df[, outgroups, drop = FALSE], get_allele))
  
  # Count how many outgroups support allele 0 or 1
  count0 <- rowSums(og == "0", na.rm = TRUE)
  count1 <- rowSums(og == "1", na.rm = TRUE)
  
  # Ancestral allele = allele supported by at least `min_support` outgroups
  ancestral <- ifelse(count0 >= min_support, "0",
                      ifelse(count1 >= min_support, "1", NA))
  
  df$ancestral <- ancestral
  
  # Keep flag: ancestral is confidently 0 (i.e. REF = ancestral allele)
  df$keep_refanc0 <- !is.na(df$ancestral) & df$ancestral == "0"
  df
}

#Select a list of outgroups and annotate ancestral state----
##outgroups <- c("AlabamaCoyote","MW692","CaliforniaCoyote") #two coyote, goldenjackal
outgroups <- c("AlabamaCoyote","MW692") #1 coyote, 1goldenjackal

total_file_anc <- add_ref_ancestral_filter(total_file, outgroups, min_support = 1) # only one supporting

# Filter to sites where ancestral is confidently 0 (REF = ancestral)
total_file_anc_real <- total_file_anc[total_file_anc$keep_refanc0, ]
# Count how many sites we lose
length(total_file_anc_real$Chrom)/length(total_file$Chrom)


#check other outgroups, not used for the final analysis, just as a sanity check
outgroups2 <- c("AlabamaCoyote","MW692","AndeanFox","BerlinZoo") #more conservative
total_file_anc2 <- add_ref_ancestral_filter(total_file, outgroups2, min_support = 3)

# Filter to sites where ancestral is confidently 0 (REF = ancestral)
total_file_anc_real2 <- total_file_anc2[total_file_anc2$keep_refanc0, ]
# Count how many sites I lose
length(total_file_anc_real2$Chrom)/length(total_file$Chrom)

# extract info for all samples (not outgroups) -----
load_list<-list()

load_list <- lapply(1:667, function(i) {
  j <- i+4
  print(i)
  
  #3 outgroups
  id <- c(samples_onlyWolves[i],outgroups,"Type")
  df <- total_file_anc_real[id]
  df$chrompos <- total_file_anc_real$chrompos

  #5 outgrouops
  #id <- c(samples_onlyWolves[i],outgroups2,"Type")
  #df <- total_file_anc_real2[id]
  #df$chrompos <- total_file_anc_real2$chrompos
  
  #metadata
  sample_name <- samples_onlyWolves[i]
  Group <- metadata_mw[metadata_mw$SampleID==samples_onlyWolves[i],]$Group
  Location <- metadata_mw[metadata_mw$SampleID==samples_onlyWolves[i],]$Location
  
  
  df_poly <- df[df[[sample_name]] != "0|0", ]
  
   colnames(df_poly) <- c("Genotype","AlabamaCoyote","MW692",
                              "Class","Variant")
  # colnames(df_poly) <- c("Genotype","AlabamaCoyote","MW692",
  #                        "AndeanFox","BerlinZoo","Class","Variant") # this for when we use more outgroups
  
   cbind.data.frame(df_poly[,c(1,4:5)], SampleID=sample_name, 
                    Group=Group, Location=Location)
  # cbind.data.frame(df_poly[,c(1,6:7)], SampleID=sample_name,
  #                  Group=Group, Location=Location)
})

library(data.table)

#load_list_genotypes <- rbindlist(load_list, use.names = TRUE, fill = TRUE)

load_list_genotypes <- do.call(rbind,load_list)

# write output  
#write.csv2(load_list_genotypes, "Files/gen_load_genotypes_wolves_only3out_NEW.csv", quote = F, row.names = F)
#write.csv2(load_list_genotypes, "Files/gen_load_genotypes_wolves_only5out_NEW.csv", quote = F, row.names = F)


# LOAD file with counts to avoid running the previous code -----
total_genotypes <- load_list_genotypes
total_genotypes <- read.csv2("Files/gen_load_genotypes_wolves_only3out_NEW.csv")
#total_genotypes <- read.csv2("Files/gen_load_genotypes_wolves_only5out_NEW.csv")

# Calculate frequencies of each variant -----

# extract derived per class 
derived_genotypes_high <- total_genotypes[total_genotypes$Genotype!="0|0"&
                                            total_genotypes$Class=="High",]
unique_variants_high <- unique(derived_genotypes_high$Variant)

derived_genotypes_moderate <- total_genotypes[total_genotypes$Genotype!="0|0"&
                                                total_genotypes$Class=="Moderate",]
unique_variants_moderate <- unique(derived_genotypes_moderate$Variant)

derived_genotypes_low <- total_genotypes[total_genotypes$Genotype!="0|0"&
                                           total_genotypes$Class=="Low",]
unique_variants_low <- unique(derived_genotypes_low$Variant)

#extract frequencies
dosage_map <- c("0|0" = 0, "0|1" = 1, "1|0" = 1, "1|1" = 2)

# High 
# n per group 
n_by_group <- derived_genotypes_high %>%
  distinct(Group, SampleID) %>%
  count(Group, name = "n_sample") %>%
  mutate(n = 2L * n_sample) %>%
  select(Group, n)

# the universe of variants to enforce
all_variants_high <- derived_genotypes_high %>%
  distinct(Variant) %>%
  arrange(Variant)

af_high <- derived_genotypes_high %>%
  mutate(dosage = unname(dosage_map[Genotype])) %>%
  group_by(Group, Variant) %>%
  summarise(sum_dosage = sum(dosage, na.rm = TRUE), .groups = "drop") %>%
  # ensure every Group has every Variant; fill missing with 0 dosage
  group_by(Group) %>%
  complete(Variant = all_variants_high$Variant, fill = list(sum_dosage = 0L)) %>%
  ungroup() %>%
  left_join(n_by_group, by = "Group") %>%
  mutate(Freq = ifelse(n > 0, sum_dosage / n, NA_real_)) %>%
  transmute(Variant, Freq, Type = "High", Group)

min_n <- 16
n_rep <- 10

sample_info <- derived_genotypes_high %>%
  filter(Group != "Zoo") %>%
  distinct(Group, SampleID)

af_down_high <- map_dfr(1:n_rep, function(rep_i) {
  
  sampled_ids <- sample_info %>%
    group_by(Group) %>%
    slice_sample(n = min_n, replace = FALSE) %>%
    ungroup()
  
  sampled_ids %>%
    crossing(Variant = all_variants_high$Variant) %>%
    left_join(
      derived_genotypes_high %>%
        filter(Group != "Zoo") %>%
        mutate(dosage = unname(dosage_map[Genotype])) %>%
        select(Group, SampleID, Variant, dosage),
      by = c("Group", "SampleID", "Variant")
    ) %>%
    mutate(dosage = replace_na(dosage, 0)) %>%
    group_by(Group, Variant) %>%
    summarise(sum_dosage = sum(dosage), .groups = "drop") %>%
    mutate(
      rep = rep_i,
      Down_freq = sum_dosage / (2 * min_n)
    )
})

dt <- as.data.table(af_down_high)
af_down_summary_high <- dt[,.(Down_freq = mean(Down_freq, na.rm = TRUE),
    SD_freq = sd(Down_freq, na.rm = TRUE),
    n_reps = sum(!is.na(Down_freq))),
  by = .(Group, Variant)]

af_high2 <- af_high %>%
  left_join(af_down_summary_high, by = c("Group", "Variant")) %>%
  mutate(Type = "High") %>%
  select(Variant, Group, Type, Freq, Down_freq, SD_freq)

# Moderate 
# n per group 
all_variants_moderate <- derived_genotypes_moderate %>%
  distinct(Variant) %>%
  arrange(Variant)

af_mod <- derived_genotypes_moderate %>%
  mutate(dosage = unname(dosage_map[Genotype])) %>%
  group_by(Group, Variant) %>%
  summarise(sum_dosage = sum(dosage, na.rm = TRUE), .groups = "drop") %>%
  # ensure every Group has every Variant; fill missing with 0 dosage
  group_by(Group) %>%
  complete(Variant = all_variants_moderate$Variant, fill = list(sum_dosage = 0L)) %>%
  ungroup() %>%
  left_join(n_by_group, by = "Group") %>%
  mutate(Freq = ifelse(n > 0, sum_dosage / n, NA_real_)) %>%
  transmute(Variant, Freq, Type = "Moderate", Group)

af_down_mod <- map_dfr(1:n_rep, function(rep_i) {
  
  sampled_ids <- sample_info %>%
    group_by(Group) %>%
    slice_sample(n = min_n, replace = FALSE) %>%
    ungroup()
  
  sampled_ids %>%
    crossing(Variant = all_variants_moderate$Variant) %>%
    left_join(
      derived_genotypes_moderate %>%
        filter(Group != "Zoo") %>%
        mutate(dosage = unname(dosage_map[Genotype])) %>%
        select(Group, SampleID, Variant, dosage),
      by = c("Group", "SampleID", "Variant")
    ) %>%
    mutate(dosage = replace_na(dosage, 0)) %>%
    group_by(Group, Variant) %>%
    summarise(sum_dosage = sum(dosage), .groups = "drop") %>%
    mutate(
      rep = rep_i,
      Down_freq = sum_dosage / (2 * min_n)
    )
})

dt <- as.data.table(af_down_mod)
af_down_summary_mod <- dt[,.(Down_freq = mean(Down_freq, na.rm = TRUE),
                              SD_freq = sd(Down_freq, na.rm = TRUE),
                              n_reps = sum(!is.na(Down_freq))),
                           by = .(Group, Variant)]

af_mod2 <- af_mod %>%
  left_join(af_down_summary_mod, by = c("Group", "Variant")) %>%
  mutate(Type = "Moderate") %>%
  select(Variant, Group, Type, Freq, Down_freq, SD_freq)


# Low 
# n per group 
all_variants_low <- derived_genotypes_low %>%
  distinct(Variant) %>%
  arrange(Variant)

af_low <- derived_genotypes_low %>%
  mutate(dosage = unname(dosage_map[Genotype])) %>%
  group_by(Group, Variant) %>%
  summarise(sum_dosage = sum(dosage, na.rm = TRUE), .groups = "drop") %>%
  # ensure every Group has every Variant; fill missing with 0 dosage
  group_by(Group) %>%
  complete(Variant = all_variants_low$Variant, fill = list(sum_dosage = 0L)) %>%
  ungroup() %>%
  left_join(n_by_group, by = "Group") %>%
  mutate(Freq = ifelse(n > 0, sum_dosage / n, NA_real_)) %>%
  transmute(Variant, Freq, Type = "Low", Group)


af_down_low <- map_dfr(1:n_rep, function(rep_i) {
  
  sampled_ids <- sample_info %>%
    group_by(Group) %>%
    slice_sample(n = min_n, replace = FALSE) %>%
    ungroup()
  
  sampled_ids %>%
    crossing(Variant = all_variants_low$Variant) %>%
    left_join(
      derived_genotypes_low %>%
        filter(Group != "Zoo") %>%
        mutate(dosage = unname(dosage_map[Genotype])) %>%
        select(Group, SampleID, Variant, dosage),
      by = c("Group", "SampleID", "Variant")
    ) %>%
    mutate(dosage = replace_na(dosage, 0)) %>%
    group_by(Group, Variant) %>%
    summarise(sum_dosage = sum(dosage), .groups = "drop") %>%
    mutate(
      rep = rep_i,
      Down_freq = sum_dosage / (2 * min_n)
    )
})

dt <- as.data.table(af_down_low)
af_down_summary_low <- dt[,.(Down_freq = mean(Down_freq, na.rm = TRUE),
                             SD_freq = sd(Down_freq, na.rm = TRUE),
                             n_reps = sum(!is.na(Down_freq))),
                          by = .(Group, Variant)]

af_low2 <- af_low %>%
  left_join(af_down_summary_low, by = c("Group", "Variant")) %>%
  mutate(Type = "Low") %>%
  select(Variant, Group, Type, Freq, Down_freq, SD_freq)



# Save frequencies
#write.csv(af_high, "Files/Frequency_derived_high_3OUT_NEW.csv", quote=F, row.names = F)
#write.csv(af_mod, "Files/Frequency_derived_moderate_3OUT_NEW.csv", quote=F, row.names = F)
#write.csv(af_low, "Files/Frequency_derived_low_3OUT_NEW.csv", quote=F, row.names = F)

write.csv(af_high2, "Files/Frequency_derived_high_3OUT_NEW_down.csv", quote=F, row.names = F)
write.csv(af_mod2, "Files/Frequency_derived_moderate_3OUT_down_NEW.csv", quote=F, row.names = F)
write.csv(af_low2, "Files/Frequency_derived_low_3OUT_NEW_down.csv", quote=F, row.names = F)

#write.csv(af_high, "Files/Frequency_derived_high_5OUT_NEW.csv", quote=F, row.names = F)
#write.csv(af_mod, "Files/Frequency_derived_moderate_5OUT_NEW.csv", quote=F, row.names = F)
#write.csv(af_low, "Files/Frequency_derived_low_5OUT_NEW.csv", quote=F, row.names = F)

# READ FREQUENCIES -----
#listfreq_moderate_df<- read.table("Files/Frequency_derived_moderate_3OUT_NEW.csv", header = T, sep=",")
#listfreq_high_df<- read.table("Files/Frequency_derived_high_3OUT_NEW.csv", header = T, sep=",")
#listfreq_low_df<- read.table("Files/Frequency_derived_low_3OUT_NEW.csv", header = T, sep=",")
#listfreq_high_df<- read.table("Files/Frequency_derived_high_5OUT_NEW.csv", header = T, sep=",")
#listfreq_moderate_df<- read.table("Files/Frequency_derived_moderate_5OUT_NEW.csv", header = T, sep=",")
#listfreq_low_df<- read.table("Files/Frequency_derived_low_5OUT_NEW.csv", header = T, sep=",")

listfreq_high_df<-af_high2
listfreq_moderate_df<-af_mod2
listfreq_low_df <- af_low2

frequency_all <- rbind.data.frame(listfreq_high_df,listfreq_moderate_df,listfreq_low_df)
frequency_all$Type <- factor(frequency_all$Type, levels=c("Low","Moderate","High"), ordered = T)
#
# RXY ASIAN -----
types <- c("High","Moderate", "Low")


frequency_all_spread_g <- frequency_all %>%
  group_by(Group, Variant, Type) %>%
  summarise(Freq = mean(Freq, na.rm = TRUE), .groups = "drop") %>%
  filter(Group != "Eastern Asian") %>%
  left_join(
    frequency_all %>%
      group_by(Group, Variant, Type) %>%
      summarise(Freq_EAS = mean(Freq, na.rm = TRUE), .groups = "drop") %>%
      filter(Group == "Eastern Asian") %>%
      select(-Group),
    by = c("Variant", "Type")
  )

Rxy_frame<-data.frame(matrix(ncol = 4, nrow = 0))
Rxy_jn<-data.frame(matrix(ncol = 3, nrow = 0))

y<-frequency_all_spread_g$Freq_EAS*(1-frequency_all_spread_g$Freq)
x<-frequency_all_spread_g$Freq*(1-frequency_all_spread_g$Freq_EAS)

xdata_total<-cbind.data.frame(x,y, Type=frequency_all_spread_g$Type, Pop=frequency_all_spread_g$Group)
Rxy_framelist <- list()

for (pop in c(groups$Group[1:12],groups$Group[14:15])){
  Rxy_frame<- data.frame()
  print(pop)
  for (j in types){
    print(j)
    xdata <- xdata_total[xdata_total$Type==j & xdata_total$Pop==pop,]
    
    # number of jackknife blocks and percentage of sites per block depend on type
    if (j == "High") {
      n_blocks <- 10
      frac <- 0.10    # 10% of the sites per block
    } else {
      n_blocks <- 100
      frac <- 0.01    # 1% of the sites per block
    }
    
    size <- floor(length(xdata$x) * frac) # take frac of the sites for the jackknife
    
    rxy <- function(a, xdata, b){
      sum(xdata[seq(a, a+(size-1), 1), b])
    } 
    
    rx_ry <- data.frame()
    for(i in c(1:n_blocks)) {
      rx_ry[i,1] <- rxy(size*(i-1)+1, xdata, 1)
      rx_ry[i,2] <- rxy(size*(i-1)+1, xdata, 2)
    } 
    
    rat <- function(c, rx_ry){ sum(rx_ry[c,1]) / sum(rx_ry[c,2]) }
    
    jn_h <- jackknife(1:n_blocks, rat, rx_ry) # do jackknife for dataframe rx_ry
    
    Rxy_jn <- rbind(Rxy_jn, cbind(jn_h$jack.values, rep(j, n_blocks)))
    Rxy_frame <- rbind(Rxy_frame, cbind(j, sum(xdata$x)/sum(xdata$y), jn_h$jack.se))
  } 
  Rxy_framelist[[pop]] <- cbind(Rxy_frame, Pop=pop)
}

Rxy_frame_df <- do.call(rbind, Rxy_framelist)
colnames(Rxy_frame_df)<-c("Type","Rxy", "SE","POP")
Rxy_frame_df$Rxy <- as.numeric(Rxy_frame_df$Rxy)
Rxy_frame_df$SE <- as.numeric(Rxy_frame_df$SE)

# Normalize by low 
Rxy_frame_low <- Rxy_frame_df[Rxy_frame_df$Type=="Low",]
Rxy_frame_notlow <- Rxy_frame_df[Rxy_frame_df$Type!="Low",]

Rxy_frame_normalized<- data.frame(Type=c("High","Moderate"),
                                  Pop=rep(c(groups$Group[1:12],groups$Group[14:15]), each=2),
                                  Rxy_norm=as.numeric(Rxy_frame_notlow$Rxy)/as.numeric(Rxy_frame_low$Rxy),
                                  SE_norm=as.numeric(Rxy_frame_notlow$S)/as.numeric(Rxy_frame_low$Rxy))

cut.off = 1
Rxy_frame_normalized$Rxy_norm_2 <- Rxy_frame_normalized$Rxy_norm - 1

Rxy_frame_normalized$Group <- Rxy_frame_normalized$Pop
Rxy_frame_normalized_m <- merge(Rxy_frame_normalized, groups[,c("Group","SuperGroup","colour","GroupID")], by="Group")

rxy_plot_all <- ggplot(Rxy_frame_normalized, aes(Type, Rxy_norm_2, fill=Type))
rxy_plot_all <- rxy_plot_all + geom_col(position = "dodge") + theme_classic()+ facet_grid(.~Pop)+
  geom_errorbar(aes(ymin=(Rxy_norm-1)-2*SE_norm, ymax=(Rxy_norm-1)+2*SE_norm), width=.2, position=position_dodge(.75)) +
  geom_hline(linetype="dashed",yintercept = 0) +scale_fill_manual(values=c("#696047","#7dde92")) +
  scale_y_continuous(labels = function(x) x + cut.off) +
  xlab("") + ylab("Rxy (Pop/Eastern Asian)") + theme(legend.position = "none",
                                                     strip.text.x = element_text(angle = 90))
rxy_plot_all

#no Asian
Rxy_frame_normalized_m$Group <- factor(Rxy_frame_normalized_m$Group, levels = order_pop2, ordered = TRUE)

col_vec <- setNames(groups_noAsian$colour, groups_noAsian$GroupID)
xlab_vec <- setNames(
  paste0("<span style='color:", col_vec, "'>", groups_noAsian$GroupID, "</span>"),
  groups_noAsian$GroupID
)

rxy_plot_noasia <- ggplot(Rxy_frame_normalized_m[-which(Rxy_frame_normalized_m$SuperGroup%in%c("Asian","Zoo")),], aes(GroupID, Rxy_norm_2, fill=Type))
rxy_plot_noasia <- rxy_plot_noasia + geom_col(position = "dodge") + theme_classic()+ facet_grid(.~SuperGroup, scales="free",space="free")+
  geom_errorbar(aes(ymin=(Rxy_norm-1)-2*SE_norm, ymax=(Rxy_norm-1)+2*SE_norm), width=.2, position=position_dodge(.75)) +
  geom_hline(linetype="dashed",yintercept = 0) + scale_fill_manual(values=c("#696047","#7dde92")) +
  #scale_fill_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  scale_y_continuous(labels = function(x) x + cut.off) +
  xlab("") + ylab("Rxy (Pop/Eastern Asian)") + 
  theme(legend.position = "bottom",
      axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5))+
  scale_x_discrete(labels = xlab_vec)
rxy_plot_noasia

ggsave("Plots/Rxy_geneticload_noAsian_3OUT_NEW.pdf", height = 5, width =12)
#ggsave("Plots/Rxy_geneticload_noAsian_5OUT_NEW.pdf", height = 5, width =12)
ggsave("Plots/FigS23.pdf", height = 5, width =12)


# RXY Scandinavian  -----
types <- c("High","Moderate", "Low")

frequency_all_spread_g <- frequency_all %>%
  group_by(Group, Variant, Type) %>%
  summarise(Freq = mean(Freq, na.rm = TRUE), .groups = "drop") %>%
  filter(Group != "Scandinavian") %>%
  left_join(
    frequency_all %>%
      group_by(Group, Variant, Type) %>%
      summarise(Freq_SCAN = mean(Freq, na.rm = TRUE), .groups = "drop") %>%
      filter(Group == "Scandinavian") %>%
      select(-Group),
    by = c("Variant", "Type")
  )

Rxy_frame<-data.frame(matrix(ncol = 4, nrow = 0))
Rxy_jn<-data.frame(matrix(ncol = 3, nrow = 0))

y<-frequency_all_spread_g$Freq_SCAN*(1-frequency_all_spread_g$Freq)
x<-frequency_all_spread_g$Freq*(1-frequency_all_spread_g$Freq_SCAN)

xdata_total<-cbind.data.frame(x,y, Type=frequency_all_spread_g$Type, Pop=frequency_all_spread_g$Group)
Rxy_framelist <- list()

for (pop in c(groups$Group[1:3],groups$Group[5:15])){
  Rxy_frame<- data.frame()
  for (j in types){
    print(j)
    xdata <- xdata_total[xdata_total$Type==j & xdata_total$Pop==pop,]
    
    # number of jackknife blocks and percentage of sites per block depend on type
    if (j == "High") {
      n_blocks <- 10
      frac <- 0.10    # 10% of the sites per block
    } else {
      n_blocks <- 100
      frac <- 0.01    # 1% of the sites per block
    }
    
    size <- floor(length(xdata$x) * frac) # take frac of the sites for the jackknife
    
    rxy <- function(a, xdata, b){
      sum(xdata[seq(a, a+(size-1), 1), b])
    } 
    
    rx_ry <- data.frame()
    for(i in c(1:n_blocks)) {
      rx_ry[i,1] <- rxy(size*(i-1)+1, xdata, 1)
      rx_ry[i,2] <- rxy(size*(i-1)+1, xdata, 2)
    } 
    
    rat <- function(c, rx_ry){ sum(rx_ry[c,1]) / sum(rx_ry[c,2]) }
    
    jn_h <- jackknife(1:n_blocks, rat, rx_ry) # do jackknife for dataframe rx_ry
    
    Rxy_jn <- rbind(Rxy_jn, cbind(jn_h$jack.values, rep(j, n_blocks)))
    Rxy_frame <- rbind(Rxy_frame, cbind(j, sum(xdata$x)/sum(xdata$y), jn_h$jack.se))
  } 
  Rxy_framelist[[pop]] <- cbind(Rxy_frame, Pop=pop)
}

Rxy_frame_df <- do.call(rbind, Rxy_framelist)


colnames(Rxy_frame_df)<-c("Type","Rxy", "SE","POP")
Rxy_frame_df$Rxy <- as.numeric(Rxy_frame_df$Rxy)
Rxy_frame_df$SE <- as.numeric(Rxy_frame_df$SE)

# Normalize by low 
Rxy_frame_low <- Rxy_frame_df[Rxy_frame_df$Type=="Low",]
Rxy_frame_notlow <- Rxy_frame_df[Rxy_frame_df$Type!="Low",]

Rxy_frame_normalized<- data.frame(Type=c("High","Moderate"),
                                  Pop=rep(c(groups$Group[1:3],groups$Group[5:15]), each=2),
                                  Rxy_norm=as.numeric(Rxy_frame_notlow$Rxy)/as.numeric(Rxy_frame_low$Rxy),
                                  SE_norm=as.numeric(Rxy_frame_notlow$S)/as.numeric(Rxy_frame_low$Rxy))

cut.off = 1
Rxy_frame_normalized$Rxy_norm_2 <- Rxy_frame_normalized$Rxy_norm - 1

Rxy_frame_normalized$Group <- Rxy_frame_normalized$Pop
Rxy_frame_normalized_m <- merge(Rxy_frame_normalized, groups[,c("Group","SuperGroup","colour","GroupID")], by="Group")

rxy_plot_all <- ggplot(Rxy_frame_normalized, aes(Type, Rxy_norm_2, fill=Type))
rxy_plot_all <- rxy_plot_all + geom_col(position = "dodge") + theme_classic()+ facet_grid(.~Pop)+
  geom_errorbar(aes(ymin=(Rxy_norm-1)-2*SE_norm, ymax=(Rxy_norm-1)+2*SE_norm), width=.2, position=position_dodge(.75)) +
  geom_hline(linetype="dashed",yintercept = 0) +scale_fill_manual(values=c("#696047","#7dde92")) +
  scale_y_continuous(labels = function(x) x + cut.off) +
  xlab("") + 
 # ylab("Rxy (Pop/Eastern Asian)") +
  ylab("Rxy (Pop/Scandinavian)") +
  theme(legend.position = "none",
                                                     strip.text.x = element_text(angle = 90))
rxy_plot_all

ggsave("Plots/Rxy_geneticload_3OUT_VSSCAN.pdf", height = 5, width =16)


#no Asian
Rxy_frame_normalized_m$Group <- factor(Rxy_frame_normalized_m$Group, levels = order_pop2, ordered = TRUE)

rxy_plot_noasia <- ggplot(Rxy_frame_normalized_m[Rxy_frame_normalized_m$SuperGroup!="Asian",], aes(Pop, Rxy_norm_2, fill=Type))
rxy_plot_noasia <- rxy_plot_noasia + geom_col(position = "dodge") + theme_classic()+ facet_grid(.~SuperGroup, scales="free",space="free")+
  geom_errorbar(aes(ymin=(Rxy_norm-1)-2*SE_norm, ymax=(Rxy_norm-1)+2*SE_norm), width=.2, position=position_dodge(.75)) +
  geom_hline(linetype="dashed",yintercept = 0) + scale_fill_manual(values=c("#696047","#7dde92")) +
  #scale_fill_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  scale_y_continuous(labels = function(x) x + cut.off) +
  xlab("") +
  ylab("Rxy (Pop/Scandinavian)") + 
  theme(legend.position = "none",axis.text.x = element_text(angle=45, hjust=1))
rxy_plot_noasia

ggsave("Plots/Rxy_geneticload_noAsian_3OUT_SCAN.pdf", height = 5, width =12)

Rxy_frame_df$Group <- Rxy_frame_df$POP
Rxy_frame_df$Rxy2 <- Rxy_frame_df$Rxy - 1
Rxy_frame_df_m <- merge(Rxy_frame_df, groups[,c("Group","SuperGroup","colour","GroupID")], by="Group")
Rxy_frame_df_m$Group <- factor(Rxy_frame_df_m$Group, levels = order_pop2, ordered = TRUE)

rxy_plot_noasia <- ggplot(Rxy_frame_df_m[Rxy_frame_df_m$SuperGroup!="Asian",], aes(Group, Rxy2, fill=Type))
rxy_plot_noasia <- rxy_plot_noasia + geom_col(position = "dodge") + theme_classic()+ facet_grid(.~SuperGroup, scales="free",space="free")+
  geom_errorbar(aes(ymin=(Rxy-1)-2*SE, ymax=(Rxy-1)+2*SE), width=.2, position=position_dodge(.75)) +
  geom_hline(linetype="dashed",yintercept = 0) + scale_fill_manual(values=c("#696047","#7dde92","grey")) +
  #scale_fill_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  scale_y_continuous(labels = function(x) x + cut.off) +
  xlab("") + ylab("Rxy (Pop/Eastern Asian)") + theme(legend.position = "none",
                                                     axis.text.x = element_text(angle=45, hjust=1))
rxy_plot_noasia


###
# Fixed variants and distributions----
df_fixed <- frequency_all %>%
  group_by(Group, Type) %>%
  summarise(
    n_total2 = n(),
    n_total  = sum(Freq > 0 , na.rm = TRUE),
    n_total_segregat  = sum(Freq > 0 & Freq<1, na.rm = TRUE),
    n_total_segregat_down  = sum(Down_freq > 0 & Down_freq<1, na.rm = TRUE),
    n_fixed_derived   = sum(Freq == 1, na.rm = TRUE),
    n_fixed_almost= sum(Freq >= 0.97, na.rm = TRUE),
    n_fixed_almost_down = sum(Down_freq >= 0.97, na.rm = TRUE),
    n_fixed_down = sum(Down_freq == 1, na.rm = TRUE),
    perc_fixed_derived   = n_fixed_derived / n_total * 100,
    perc_fixed_almostdown_derived = n_fixed_almost_down / n_total * 100
  )

df_fixed_m <- merge(df_fixed, groups[,c("Group","colour","SuperGroup","GroupID")], by="Group")
df_fixed_m$SuperGroup <- factor(df_fixed_m$SuperGroup, levels=c("Eastern","Fennoscandia","Southern","Asian","Zoo"), ordered = T)
df_fixed_m$Type <- factor(df_fixed_m$Type, levels=c("High","Moderate","Low"), ordered = T)
df_fixed_m_noasian <- df_fixed_m[-which(df_fixed_m$SuperGroup%in%c("Asian","Zoo")),]

col_vec <- setNames(groups_noAsian$colour, groups_noAsian$GroupID)
xlab_vec <- setNames(
  paste0("<span style='color:", col_vec, "'>", groups_noAsian$GroupID, "</span>"),
  groups_noAsian$GroupID
)


g_fixed <- ggplot(df_fixed_m_noasian[df_fixed_m_noasian$Type!="Low",], aes(GroupID,n_fixed_derived, fill=GroupID))
g_fixed <- g_fixed +geom_col() + theme_classic() + facet_grid(Type~SuperGroup, scales="free") +
 # geom_point(data=df_fixed_m_noasian[df_fixed_m_noasian$Type!="Low",], aes(GroupID,n_fixed_down, fill=GroupID))+
  scale_fill_manual(values = col_vec) +
  guides(color=guide_legend(title="Population"))+ xlab("")+ylab("Fixed Deleterious Derived Variants ")+
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none", strip.text=element_text(size=8))+
  scale_x_discrete(labels = xlab_vec)

g_fixed  


g_genotypes <- ggplot(df_fixed_m_noasian[df_fixed_m_noasian$Type!="Low",], aes(GroupID,n_total_segregat, fill=GroupID))
g_genotypes <- g_genotypes +geom_col() + theme_classic() + facet_grid(Type~SuperGroup, scales="free") +
  scale_fill_manual(values = col_vec) +
  #geom_point(data=df_fixed_m_noasian[df_fixed_m_noasian$Type!="Low",], aes(GroupID,n_total_segregat_down, fill=GroupID))+
  guides(color=guide_legend(title="Population"))+ xlab("")+ylab("Segregating Deleterious Variants ")+
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none", strip.text=element_text(size=8))+
  scale_x_discrete(labels = xlab_vec)

g_genotypes  



# with asian
col_vec_as <- setNames(groups$colour, groups$GroupID)
xlab_vec_as <- setNames(
  paste0("<span style='color:", col_vec_as, "'>", groups$GroupID, "</span>"),
  groups$GroupID
)

g_fixed_a <- ggplot(df_fixed_m[df_fixed_m$Type!="Low",], aes(GroupID,n_fixed_derived, fill=GroupID))
g_fixed_a <- g_fixed_a +geom_col() + theme_classic() + facet_grid(Type~SuperGroup, scales="free", space="free_x") +
  geom_point(data=df_fixed_m[df_fixed_m$Type!="Low"&df_fixed_m$Group!="Zoo",], aes(GroupID,n_fixed_down, fill=GroupID))+
  scale_fill_manual(values = col_vec_as) +
  guides(color=guide_legend(title="Population"))+ xlab("")+ylab("Fixed Deleterious Derived Variants ")+
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none")+
  scale_x_discrete(labels = xlab_vec_as)

g_fixed_a 


g_genotypes_a <- ggplot(df_fixed_m[df_fixed_m$Type!="Low",], aes(GroupID,n_total_segregat, fill=GroupID))
g_genotypes_a <- g_genotypes_a +geom_col() + theme_classic() + facet_grid(Type~SuperGroup, scales="free", space="free_x") +
  scale_fill_manual(values = col_vec_as) +
  geom_point(data=df_fixed_m[df_fixed_m$Type!="Low"&df_fixed_m$Group!="Zoo",], aes(GroupID,n_total_segregat_down, fill=GroupID))+
  guides(color=guide_legend(title="Population"))+ xlab("")+ylab("Segregating Deleterious Variants ")+
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none")+
  scale_x_discrete(labels = xlab_vec_as)

g_genotypes_a  

g_fixed_a / g_genotypes_a + plot_layout(heights = c(1, 1)) + plot_annotation(tag_levels = 'A')

ggsave("Plots/FigS22.pdf", width = 8, height = 8)



g2 <- ggplot(df_fixed_m_noasian, aes(GroupID,n_fixed_almost, fill=GroupID))
g2 <- g2 +geom_col() + theme_classic() + facet_grid(Type~SuperGroup, scales="free") +
  scale_fill_manual(values = col_vec) +
  guides(color=guide_legend(title="Population"))+ xlab("")+ylab("Fixed Derived Alleles")+
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none")+
  scale_x_discrete(labels = xlab_vec)

g2  

# Deleterious Allele COUNTS ----
#extract counts no need to run again----
types <- c("High", "Moderate","Low")
total_genotypes <- load_list_genotypes


df_counts <- total_genotypes %>%
  # keep only wolves + relevant classes
  filter(SampleID %in% samples_onlyWolves,
         Class %in% types) %>%
  group_by(SampleID, Class) %>%
  summarise(
    Group    = first(Group),
    Location = first(Location),
    Het      = sum(Genotype %in% c("0|1", "1|0")),
    Hom      = sum(Genotype == "1|1") * 2L,
    .groups  = "drop"
  )

df_counts$Total <-df_counts$Het + df_counts$Hom

#write.csv2(df_counts,"Files/Counts_all_3OUT_NEW.csv", quote=F, row.names=F)
#write.csv2(df_counts,"Files/Counts_all_5OUT_NEW.csv", quote=F, row.names=F)

# Read saved file for COUNTs ----
#df_counts <- read.csv2("Files/Counts_all_3OUT_NEW.csv")
df_counts <- read.csv2("Files/Counts_all_5OUT_NEW.csv")

df_counts_gather <- gather(df_counts, Load, Count, -SampleID, -Class,
                           -Group, -Location)

# Plot counts
df_counts_gather$Load <- gsub("Het","Heterozygous Count",df_counts_gather$Load)
df_counts_gather$Load <- gsub("Hom","Homozygous Count",df_counts_gather$Load)
df_counts_gather$Load <- gsub("Total","Total Count",df_counts_gather$Load)
df_counts_gather$Class <- factor(df_counts_gather$Class, ordered = T, levels=c("Low","Moderate","High"))

df_counts_gather_m <- merge(df_counts_gather, groups[,c("Group","colour","SuperGroup")], by="Group")

df_counts_gather_m_noAsian <- df_counts_gather_m[-which(df_counts_gather_m$SuperGroup%in%c("Asian","Zoo")),]


g_counts4 <- ggplot(df_counts_gather_m_noAsian, 
                    aes(Load, Count, color=colour))
g_counts4 <- g_counts4 + geom_boxplot(outlier.shape = NA, position = position_dodge(width=0.6)) + 
  geom_point(position=position_jitterdodge(dodge.width = 0.6, jitter.width = 0.3), size=0.2)+theme_classic() + 
  facet_grid(Class~SuperGroup, scales="free", space="free_x") + xlab("") + 
  theme(legend.position = "top", axis.text.x = element_text(angle=45, hjust=1))+
  ylab("Deleterious Allele Counts")+
  scale_color_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  guides(color=guide_legend(title="Population"))
g_counts4

#ggsave("Plots/FigS21_5OUT.pdf", width = 8, height = 7)

g_counts4 <- ggplot(df_counts_gather_m_noAsian[df_counts_gather_m_noAsian$Load!="Total Count",], 
                    aes(Load, Count, color=colour))
g_counts4 <- g_counts4 + geom_boxplot(outlier.shape = NA, position = position_dodge(width=0.6)) + 
  geom_point(position=position_jitterdodge(dodge.width = 0.6, jitter.width = 0.3), size=0.2)+theme_classic() + 
  facet_grid(Class~SuperGroup, scales="free", space="free_x") + xlab("") + 
  theme(legend.position = "top", axis.text.x = element_text(angle=45, hjust=1))+
  ylab("Deleterious Allele Counts")+
  scale_color_identity(guide = "legend", labels = groups$Group,breaks = groups$colour) +
  guides(color=guide_legend(title="Population"))
g_counts4

ggsave("Plots/FigS24.pdf", width = 8, height = 7)


# MAIN FIGURE----------
# Directly from the metadata file
metadata_mw2_noasian <- metadata_mw2[-which(metadata_mw2$SuperGroup%in%c("Asian","Zoo")),]
metadata_mw2_noasian_spcounts <- metadata_mw2_noasian %>%
  pivot_longer(
    cols = c(Total_High, Total_Moderate,Total_Low),
    names_to  = "Class",
    values_to = "TotalCounts"
  )%>%
  mutate(
    Class = recode(Class,
                   "Total_Moderate" = "Moderate",
                   "Total_High"     = "High",
                   "Total_Low"="Low")
  )
metadata_mw2_noasian_spcounts$Class<- factor(metadata_mw2_noasian_spcounts$Class, 
                                             levels=c("High","Moderate","Low"), ordered = T)
col_vec <- setNames(groups_noAsian$colour, groups_noAsian$GroupID)
xlab_vec <- setNames(
  paste0("<span style='color:", col_vec, "'>", groups_noAsian$GroupID, "</span>"),
  groups_noAsian$GroupID
)

g_counts_total <- ggplot(metadata_mw2_noasian_spcounts[metadata_mw2_noasian_spcounts$Class !="Low",], 
                         aes(GroupID, TotalCounts, color=GroupID))
g_counts_total <- g_counts_total + 
  geom_boxplot(outlier.shape = NA, position = position_dodge(width=0.6)) + 
  geom_point(position=position_jitterdodge(dodge.width = 0.6, jitter.width = 0.3), size=0.2)+theme_classic() + 
  facet_grid(Class~SuperGroup, scales="free", space="free_x") + xlab("") + 
  theme(legend.position = "none", axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  #stat_compare_means() +
  scale_color_manual(values = col_vec) +
  scale_shape_manual(values=c(1,8,17))+
  ylab("Total Deleterious Allele Counts")+
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),, 
        legend.position = "none", strip.text=element_text(size=8))+
  scale_x_discrete(labels = xlab_vec)
g_counts_total



metadata_mw2_noasian_hom_counts <- metadata_mw2_noasian %>%
  pivot_longer(
    cols = c(Hom_High, Hom_Moderate,Hom_Low),
    names_to  = "Class",
    values_to = "HomCounts"
  )%>%
  mutate(
    Class = recode(Class,
                   "Hom_Moderate" = "Moderate",
                   "Hom_High"     = "High",
                   "Hom_Low"="Low")
  )
metadata_mw2_noasian_hom_counts$Class<- factor(metadata_mw2_noasian_hom_counts$Class, 
                                             levels=c("High","Moderate","Low"), ordered = T)

g_counts_hom <- ggplot(metadata_mw2_noasian_hom_counts[metadata_mw2_noasian_hom_counts$Class !="Low",], 
                         aes(GroupID, HomCounts, color=GroupID))
g_counts_hom <- g_counts_hom + geom_boxplot(outlier.shape = NA, position = position_dodge(width=0.6)) + 
  geom_point(position=position_jitterdodge(dodge.width = 0.6, jitter.width = 0.3), size=0.2)+theme_classic() + 
  facet_grid(Class~SuperGroup, scales="free", space="free_x") + xlab("") + 
  theme(legend.position = "none", axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  #stat_compare_means() +
  scale_color_manual(values = col_vec) +
  scale_shape_manual(values=c(1,8,17))+
  ylab("Homozygous Deleterious Allele Counts")+
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),, 
        legend.position = "none", strip.text=element_text(size=8))+
  scale_x_discrete(labels = xlab_vec)
g_counts_hom


#correlations Froh and homozygous load
corr <- ggplot(metadata_mw2_noasian_hom_counts[metadata_mw2_noasian_hom_counts$Class !="Low",], 
               aes(FROH, HomCounts, color=GroupID))
corr <- corr + 
  geom_point(alpha=0.4)+theme_classic() +
  geom_smooth(method = "lm")+
  facet_grid(Class~SuperGroup, scales="free_y") + xlab(expression(F[ROH])) + 
  theme(legend.position = "none", axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  scale_color_manual(values = col_vec) +
  scale_shape_manual(values=c(1,8,17))+
  ylab("Homozygous Deleterious Allele Counts")
corr


(g_counts_total + g_genotypes ) / (g_counts_hom + g_fixed)  + plot_annotation(tag_levels = 'A')
ggsave("Plots/Figure4_new.pdf", height = 17, width = 18.4,units="cm")

# Some stats -----
metadata_mw2 %>%
  group_by(Group) %>%
  summarise(
    mean_FROH = mean(FROH, na.rm = TRUE),
    sd_FROH   = sd(FROH, na.rm = TRUE),        # optional
    n_samples = n()                             # optional
  ) %>% arrange(desc(mean_FROH))

metadata_mw2 %>%
  group_by(Group) %>%
  summarise(
    mean_Heterozygosity = mean(Het_VCF, na.rm = TRUE),
    sd_Heterozygosity   = sd(Het_VCF, na.rm = TRUE),        # optional
    n_samples = n()                             # optional
  ) %>% arrange(desc(mean_Heterozygosity))


# increase in load - maybe remove Zoo as well
load_summ <- metadata_mw2 %>%
  filter(!Group %in% c("Eastern Asian","South Western Asian")) %>%
  group_by(Group) %>%
  summarise(
    ratio_Total_Moderate = mean(Total_Moderate, na.rm = TRUE)/mean(metadata_mw2[which(!metadata_mw2$Group %in% c("Eastern Asian", "South Western Asian")),]$Total_Moderate),
    perc_Total_Moderate = (mean(Total_Moderate, na.rm = TRUE)/mean(metadata_mw2[which(!metadata_mw2$Group %in% c("Eastern Asian","South Western Asian")),]$Total_Moderate)-1)*100,
    mean_Total_Moderate = mean(Total_Moderate, na.rm = TRUE),
    
    ratio_Total_High = mean(Total_High, na.rm = TRUE)/mean(metadata_mw2[which(!metadata_mw2$Group %in% c("Eastern Asian","South Western Asian")),]$Total_High),
    perc_Total_High = (mean(Total_High, na.rm = TRUE)/mean(metadata_mw2[which(!metadata_mw2$Group %in% c("Eastern Asian", "South Western Asian")),]$Total_High)-1)*100,
    mean_Total_High = mean(Total_High, na.rm = TRUE),
    
    #ratio_T_Moderate = Total_Moderate/mean_Total_Moderate,
    #ratio_T_High = Total_High/mean_Total_High,
    n_samples = n()                             # optional
  )# %>% arrange(desc(mean_Total_High))
load_summ
write.csv2(load_summ,"Files/Summary_geneticload_average.csv", quote = F, row.names = F)

# segregating mutation
sites_high <- df_fixed_m %>%
  filter(!Group %in% c("Eastern Asian","South Western Asian")) %>%
  filter(Type %in% c("High")) %>%
  group_by(Group) %>%
  reframe(
    ratio_Total_Segregating = n_total_segregat/mean(df_fixed_m[which(!df_fixed_m$Group %in% c("Eastern Asian","South Western Asian")&df_fixed_m$Type %in% c("High")),]$n_total_segregat), 
    perc_Total_Segregating = (n_total_segregat/mean(df_fixed_m[which(!df_fixed_m$Group %in% c("Eastern Asian","South Western Asian")&df_fixed_m$Type %in% c("High")),]$n_total_segregat)-1)*100,
    mean_Total_Segregating = n_total_segregat,
    
    ratio_Total_fixed_derived = n_fixed_derived/mean(df_fixed_m[which(!df_fixed_m$Group %in% c("Eastern Asian","South Western Asian")&df_fixed_m$Type %in% c("High")),]$n_fixed_derived),
    perc_Total_fixed_derived = (n_fixed_derived/mean(df_fixed_m[which(!df_fixed_m$Group %in% c("Eastern Asian","South Western Asian")&df_fixed_m$Type %in% c("High")),]$n_fixed_derived)-1)*100,
    mean_Total_fixed_derived = n_fixed_derived,
    
    n_samples = n()                             # optional
  ) %>% arrange(ratio_Total_Segregating)
sites_high
write.csv2(sites_high,"Files/Summary_sites_high_average.csv", quote = F, row.names = F)


sites_mod <- df_fixed_m %>%
  filter(!Group %in% c("Eastern Asian", "South Western Asian")) %>%
  filter(Type %in% c("Moderate")) %>%
  group_by(Group) %>%
  reframe(
    ratio_Total_Segregating = n_total_segregat/mean(df_fixed_m[which(!df_fixed_m$Group %in% c("Eastern Asian", "Volga-Ural","South Western Asian")&df_fixed_m$Type %in% c("Moderate")),]$n_total_segregat), 
    perc_Total_Segregating = (n_total_segregat/mean(df_fixed_m[which(!df_fixed_m$Group %in% c("Eastern Asian", "Volga-Ural","South Western Asian")&df_fixed_m$Type %in% c("Moderate")),]$n_total_segregat)-1)*100,
    mean_Total_Segregating = n_total_segregat,
    
    ratio_Total_fixed_derived = n_fixed_derived/mean(df_fixed_m[which(!df_fixed_m$Group %in% c("Eastern Asian", "Volga-Ural","South Western Asian")),]$n_fixed_derived),
    perc_Total_fixed_derived = (n_fixed_derived/mean(df_fixed_m[which(!df_fixed_m$Group %in% c("Eastern Asian", "Volga-Ural","South Western Asian")),]$n_fixed_derived)-1)*100,
    mean_Total_fixed_derived = n_fixed_derived,

    n_samples = n()                             # optional
  ) %>% arrange(ratio_Total_Segregating)

write.csv2(sites_mod,"Files/Summary_sites_moderate_average.csv", quote = F, row.names = F)
