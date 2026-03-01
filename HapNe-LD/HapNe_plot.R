#HapNe trajectories - Recent Ne
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(readxl)
library(ggrepel)
library(ggtext)
library(patchwork)

setwd("~/Documents/OneDrive - University of Copenhagen/ResQ_MSCA/")

#Metadata -----
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


# Read files recent Ne with HapNe LD all samples------

# 1) List all csv files from HapNe_try1 ... HapNe_try10
ne_files <- unlist(lapply(1:10, function(i) {
  list.files(paste0("Files/HapNe_all/HapNe_try", i, "/"),
             pattern = "\\.csv$",
             full.names = TRUE)
}))


groups_ne <- sub('\\ld_hapne_estimate.csv$', '', basename(ne_files))
groups_ne <- gsub('M_','', groups_ne)
groups_ne <- gsub('_',' ', groups_ne)

groups_ne <- gsub('NorthWesternIberian','North Western Iberian', groups_ne)
groups_ne <- gsub('Southern Western Asian','South Western Asian', groups_ne)
groups_ne <- gsub('WesternFinland','Western Finland', groups_ne)


try_ne <- sub(".*HapNe_try([0-9]+).*", "\\1", ne_files)

ne_list <- lapply(ne_files,FUN = function(files) {
  read.csv2(files, sep=",")
})

# Bind into one df with Group + try + content
ne_list2 <- vector("list", length(ne_files))
for (i in seq_along(ne_files)) {
  ne_list2[[i]] <- cbind.data.frame(
    Group = groups_ne[i],
    Rep   = as.integer(try_ne[i]),
    ne_list[[i]]
  )
}

# Final dataframe
ne_df <- do.call(rbind, ne_list2)



#Plotting -----

ne_df$Group <- factor(ne_df$Group, levels = order_pop, ordered = T)
ne_df_m <- merge(ne_df, groups[,c("Group","colour","SuperGroup")], by="Group")


# average all rep
ne_avg <- ne_df_m %>%
  group_by(Group, TIME) %>%
  summarise(
    Q0.5_median = median(Q0.5, na.rm = TRUE),
    Q0.5_low  = quantile(Q0.5, 0.025, na.rm = TRUE),
    Q0.5_high = quantile(Q0.5, 0.975, na.rm = TRUE),
    .groups = "drop"
  )
ne_avg <- merge(ne_avg, groups[,c("Group","colour","SuperGroup")], by="Group")



g_g <- ggplot() +
  # Background: all other tries in grey
  geom_line(
    data = ne_df_m,
    aes(x = TIME, y = Q0.5, group = interaction(Group, Rep)),
    color = "grey75", linewidth = 0.7, alpha = 0.7
  ) +
  # Average
  #geom_line(data=ne_avg, aes(x = TIME, y = Q0.5_median, color = colour, group = Group), linewidth = 1) +
  geom_smooth(data = ne_avg,aes(x = TIME, y = Q0.5_median, group = Group, color = colour),method = "gam",
    formula = y ~ s(x, k = 10), se = FALSE,  linewidth = 1.1)  +
  # Confidence ribbons ONLY from average 1
  geom_ribbon(
    data = ne_avg,
    aes(x = TIME, ymin = Q0.5_low, ymax = Q0.5_high, group = Group, fill=colour),
     color = NA, alpha = 0.2
  ) +
  facet_wrap(.~Group, scales = "free_y", nrow = 3) +
  labs(x = "Generations ago", y = "Ne") +
  scale_color_identity(
    guide  = "legend",
    labels = groups$Group,
    breaks = groups$colour
  ) +
  scale_fill_identity(
    guide  = "legend",
    labels = groups$Group,
    breaks = groups$colour
  ) +
  guides(color = guide_legend(title = "Population")) +
 # coord_cartesian(clip = "on", ylim = c(0,15000)) +
  scale_x_continuous(
    limits = c(0, 50),
    n.breaks = 10,
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_classic() +
  theme(legend.position = "none")

g_g

ggsave("Plots/HapNeLd_50g_allsamples.pdf", height = 8, width = 14)

# Summary ----
summaryNe <- list()
for (i in 1:length(order_pop)){
  minNe <- min(ne_avg[ne_avg$Group==order_pop[i],]$Q0.5_median)
  maxNe <- max(ne_avg[ne_avg$Group==order_pop[i],]$Q0.5_median)
  Ne50 <- min(ne_avg[ne_avg$Group==order_pop[i]&ne_avg$TIME==50,]$Q0.5_median)
  Ne20 <- max(ne_avg[ne_avg$Group==order_pop[i]&ne_avg$TIME==20,]$Q0.5_median)
  Ne10 <- max(ne_avg[ne_avg$Group==order_pop[i]&ne_avg$TIME==10,]$Q0.5_median)
  Ne1 <- max(ne_avg[ne_avg$Group==order_pop[i]&ne_avg$TIME==1,]$Q0.5_median)
  Group <- order_pop[i]
  Supergroup <- ne_df_m[ne_df_m$Group==order_pop[i],]$SuperGroup[1]
  summaryNe[[i]] <- data.frame(Group=Group, SuperGroup=Supergroup, minNe=minNe,maxNe=maxNe,
                               Ne50=Ne50,Ne20=Ne20,Ne10=Ne10, Ne1=Ne1)
  
}

summaryNe_df <- do.call(rbind,summaryNe)


#Add rate of decline (or increase) - Exponential rate of change
summaryNe_df$r_50_20 <- (log(summaryNe_df$Ne50) - log(summaryNe_df$Ne20)) / (50 - 20)
summaryNe_df$r_20_10 <- (log(summaryNe_df$Ne20) - log(summaryNe_df$Ne10)) / (20 - 10)
summaryNe_df$r_10_1  <- (log(summaryNe_df$Ne10)  - log(summaryNe_df$Ne1)) / (10 - 1)
summaryNe_df$r_50_10  <- (log(summaryNe_df$Ne50)  - log(summaryNe_df$Ne1)) / (50 - 10)

summaryNe_df$perc_50_20 <- (exp(summaryNe_df$r_50_20) - 1) * 100
summaryNe_df$perc_20_10  <- (exp(summaryNe_df$r_20_10) - 1) * 100
summaryNe_df$perc_10_1  <- (exp(summaryNe_df$r_10_1) - 1) * 100
summaryNe_df$perc_50_10  <- (exp(summaryNe_df$r_50_10) - 1) * 100


write.csv2(summaryNe_df,"Files/Ne_rates.csv", quote=F, row.names = F)

summaryNe_df_g <- gather(summaryNe_df, Rate_G, Perc, -Group, -SuperGroup, -minNe,-maxNe,-Ne50,
                                        -Ne20,-Ne10,-Ne1,-r_50_20,-r_20_10,-r_10_1,-r_50_10)

summaryNe_df_gm <- merge(summaryNe_df_g, groups[,c("Group","colour","GroupID")], by="Group")

summaryNe_df_gm$Group <- factor(summaryNe_df_gm$Group, ordered = T, levels = order_pop)

summaryNe_df_gm$Rate_G <- gsub("perc_50_10","Gen50_10",summaryNe_df_gm$Rate_G )
summaryNe_df_gm$Rate_G <- gsub("perc_20_10","Gen20_10",summaryNe_df_gm$Rate_G )
summaryNe_df_gm$Rate_G <- gsub("perc_10_1","Gen10_1",summaryNe_df_gm$Rate_G )
summaryNe_df_gm$Rate_G <- gsub("perc_50_20","Gen50_20",summaryNe_df_gm$Rate_G )

summaryNe_df_gm$Rate_G <- factor(summaryNe_df_gm$Rate_G, levels=c("Gen10_1","Gen20_10","Gen50_20","Gen50_10"), ordered = T)
summaryNe_df_gm2 <- summaryNe_df_gm[summaryNe_df_gm$Rate_G!="Gen50_10",]


groups_noAsian <- groups[groups$Group %in% order_pop2, ]
groups_noAsian <- groups_noAsian[match(order_pop2, groups_noAsian$Group), ]
summaryNe_df_gm_noAsian <- summaryNe_df_gm2[summaryNe_df_gm2$SuperGroup!="Asian",]

col_vec <- setNames(groups_noAsian$colour, groups_noAsian$GroupID)
xlab_vec <- setNames(
  paste0("<span style='color:", col_vec, "'>", groups_noAsian$GroupID, "</span>"),
  groups_noAsian$GroupID
)

g_rate <- ggplot(summaryNe_df_gm_noAsian, aes(GroupID,Perc*-1, fill=Group)) #make it that negative is decline
g_rate <- g_rate + geom_col() + facet_grid(Rate_G~SuperGroup, scales="free",space="free_x")+
  scale_fill_manual(
    values = setNames(groups$colour, groups$Group),  # hex codes mapped to Group
    breaks = groups$Group,                           # ordering of legend
    labels = groups$Group                            # legend labels
  )+
  guides(fill=guide_legend(title="Population"))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_classic() + 
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust=0.5),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.8),
        legend.position = "none")+  scale_x_discrete(labels = xlab_vec)+
  xlab("") + ylab("Percent of exponential change per generation")
  
g_rate

ggsave("Plots/FigS19.pdf", width = 8, height = 6)


summaryNe_df_gm_one <- summaryNe_df_gm[summaryNe_df_gm$Rate_G=="Gen50_10",]
summaryNe_df_gm_one_noAsian <- summaryNe_df_gm_one[summaryNe_df_gm_one$SuperGroup!="Asian",]

g_rate2 <- ggplot(summaryNe_df_gm_one_noAsian, aes(GroupID,Perc*-1, fill=Group)) #make it that negative is decline
g_rate2 <- g_rate2 + geom_col() + facet_grid(Rate_G~SuperGroup, scales="free",space="free_x")+
  scale_fill_manual(
    values = setNames(groups$colour, groups$Group),  # hex codes mapped to Group
    breaks = groups$Group,                           # ordering of legend
    labels = groups$Group                            # legend labels
  )+
  guides(fill=guide_legend(title="Population"))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_classic() + 
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust=0.5),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.8),
        legend.position = "none")+  scale_x_discrete(labels = xlab_vec)+
  xlab("") + ylab("Percent of exponential change per generation")

g_rate2


g <- ggplot(summaryNe_df_gm, aes(Group,Perc*-1, fill=colour)) #make it that negative is decline
g <- g + geom_col() + facet_grid(.~Rate_G, scales="free",space="free_x")+
  scale_fill_identity(guide = "legend", labels = groups$Group, breaks = groups$colour) +
  guides(fill=guide_legend(title="Population"))+
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "none") +
  xlab("") + ylab("Percent of exponential change per generation")

g

# only ASIAN
lab_df <- ne_avg[ne_avg$SuperGroup =="Asian"&ne_avg$TIME==50,] 

#50 generations ago and supergorup
g_gAsian <- ggplot()
g_gAsian <- g_gAsian + 
  # Background: all tries in grey
  geom_line(data = ne_df_m[ne_df_m$SuperGroup=="Asian",], aes(x = TIME, y = Q0.5, group = interaction(Group, Rep), color=colour),
            #stat="smooth", method = "gam", formula = y ~ s(x, k = 10),
            linewidth = 0.5,
            alpha = 0.5)+
  # Average
  geom_line(data = ne_avg[ne_avg$SuperGroup=="Asian",],aes(x = TIME, y = Q0.5_median, group = Group, color = colour),
            # stat="smooth", method = "gam", formula = y ~ s(x, k = 10), se = FALSE,  
            linewidth = 1.1) +
  facet_wrap(.~SuperGroup, scales="free_y",nrow=1)+
  xlab(paste0("Generations ago")) + ylab("Ne")+
  scale_color_identity(guide = "legend", labels = groups$Group, breaks = groups$colour) +
  scale_fill_identity(guide = "legend", labels = groups$Group, breaks = groups$colour) +
  guides(color=guide_legend(title="Population"))+
  scale_x_continuous(limits=c(0,50),
                     n.breaks = 10,labels = scales::comma,expand=expansion(mult = c(0, 0.05)))+
  theme_classic()+ theme(legend.position = "none") +
  geom_text_repel(data = lab_df,
    aes(x = TIME, y = Q0.5_median,label = Group, color=colour),
    size = 3,
    nudge_x = -1,         # move slightly left from edge
    direction = "y",      # only move vertically
    box.padding = 0.3,
    segment.color = "grey70",
    segment.size = 0.3,
    max.overlaps = Inf,
    show.legend = FALSE )
g_gAsian

ggsave("Plots/HapNeLD_onlyAsian.pdf", width = 5, height = 3)
ggsave("Plots/FigS18.pdf", width = 5, height = 3)



#Make main Figure 3 ------

# PLOT by SUBGROUP without Asians
df_filt <- ne_df_m %>%
  filter(SuperGroup != "Asian", TIME <= 50)

# Only average
lab_df <- ne_avg %>%
  group_by(Group, SuperGroup) %>%
  filter(SuperGroup != "Asian", TIME <= 50)%>%
  ungroup() %>%
  mutate(
    Q0.5 = Q0.5_median   # rename so geom_text_repel works without changes
  )
# base plot
g_g0 <- ggplot()
g_g0 <- g_g0 + 
   # Background: all tries in grey
  geom_line(data = df_filt, aes(x = TIME, y = Q0.5, group = interaction(Group, Rep), color=colour),
           #stat="smooth", method = "gam", formula = y ~ s(x, k = 10),
            linewidth = 0.5,
            alpha = 0.5)+
  # Average
  geom_line(data = ne_avg[ne_avg$SuperGroup!="Asian",],aes(x = TIME, y = Q0.5_median, group = Group, color = colour),
           # stat="smooth", method = "gam", formula = y ~ s(x, k = 10), se = FALSE,  
            linewidth = 1.1)  +
  facet_wrap(.~SuperGroup, scales="free_y",nrow=1)+
  xlab(paste0("Generations ago")) + ylab("Ne")+
  scale_color_identity(guide = "legend", labels = groups$Group, breaks = groups$colour) +
  scale_fill_identity(guide = "legend", labels = groups$Group, breaks = groups$colour) +
  guides(color=guide_legend(title="Population"))+
  
  coord_cartesian(clip = "on", ylim = c(0,10000))+
  scale_x_continuous(limits=c(0,50),
                     n.breaks = 10,labels = scales::comma,expand=expansion(mult = c(0, 0.05)))+
  theme_classic()+ theme(legend.position = "none") +
  geom_text_repel(data = lab_df[lab_df$TIME==50,],
                  aes(TIME, Q0.5_median,label = Group, color=colour),
                  size = 3,
                  nudge_x = -1,         # move slightly left from edge
                  nudge_y = -1,         # move slightly left from edge
                  direction = "y",      # only move vertically
                  box.padding = 0.3,
                  segment.color = "grey70",
                  segment.size = 0.3,
                  max.overlaps = Inf,
                  show.legend = FALSE
  )
g_g0
ggsave("Plots/HapNeLd_50g_allsamples_trajectory_allruns.pdf", height = 8, width = 14)


# inserts 
g_g1 <- ggplot( )
g_g1 <- g_g1 + 
  geom_line(data = df_filt[df_filt$SuperGroup=="Eastern"&df_filt$Group!="Eastern European Lowlands",], aes(x = TIME, y = Q0.5, group = interaction(Group, Rep), color=colour),
            #stat="smooth",method = "gam", formula = y ~ s(x, k = 10),
            linewidth = 0.5,
            alpha = 0.3)+
  geom_line(data=lab_df[lab_df$SuperGroup=="Eastern"&lab_df$TIME<=10&lab_df$Group!="Eastern European Lowlands",],
            aes(x = TIME, y = Q0.5_median, color=colour),
            linewidth=1.1)+
  xlab("") + ylab("")+
  coord_cartesian(ylim = c(0, 1000)) + 
  scale_color_identity(guide = "legend", labels = groups$Group, breaks = groups$colour) +
  guides(color=guide_legend(title="Population"))+
  scale_x_continuous(limits=c(0,10),
                     n.breaks = 10,labels = scales::comma,expand=expansion(mult = c(0, 0.05)))+
  theme_classic()+ theme(legend.position = "none")
g_g1

g_g2 <- ggplot()
g_g2 <- g_g2 + 
  geom_line(data = df_filt[df_filt$SuperGroup=="Fennoscandia"&df_filt$Group!="Karelian",], aes(x = TIME, y = Q0.5, group = interaction(Group, Rep), color=colour),
           # stat="smooth",method = "gam", formula = y ~ s(x, k = 10),
            linewidth = 0.5,
            alpha = 0.3)+
  geom_line(data=lab_df[lab_df$SuperGroup=="Fennoscandia"&lab_df$TIME<=10&lab_df$Group!="Karelian",],aes(x = TIME, y = Q0.5_median, color=colour),
            linewidth=1.1)+
  xlab("") + ylab("")+
  coord_cartesian(ylim = c(0, 300)) +
  scale_color_identity(guide = "legend", labels = groups$Group, breaks = groups$colour) +
  guides(color=guide_legend(title="Population"))+
  scale_x_continuous(limits=c(0,10),
                     n.breaks = 10,labels = scales::comma,expand=expansion(mult = c(0, 0.05)))+
  theme_classic()+ theme(legend.position = "none")
g_g2

g_g3 <- ggplot(lab_df[lab_df$SuperGroup=="Southern"&lab_df$TIME<=10,],aes(x = TIME, y = Q0.5, color=colour) )
g_g3 <- g_g3 + 
  geom_line(data = df_filt[df_filt$SuperGroup=="Southern",], aes(x = TIME, y = Q0.5, group = interaction(Group, Rep), color=colour),
          #  stat="smooth",method = "gam", formula = y ~ s(x, k = 10),
            linewidth = 0.5,
            alpha = 0.3)+
  geom_line(data=lab_df[lab_df$SuperGroup=="Southern"&lab_df$TIME<=10,],aes(x = TIME, y = Q0.5_median, color=colour),
            linewidth=1.1)+
  xlab("") + ylab("")+
  coord_cartesian(ylim = c(0, 1500)) + 
  scale_color_identity(guide = "legend", labels = groups$Group, breaks = groups$colour) +
  guides(color=guide_legend(title="Population"))+
  scale_x_continuous(limits=c(0,10),
                     n.breaks = 10,labels = scales::comma,expand=expansion(mult = c(0, 0.05)))+
  theme_classic()+ theme(legend.position = "none")
g_g3


#function to annotate
annotation_custom2 <- function(grob, xmin = -Inf, xmax = Inf,
                               ymin = -Inf, ymax = Inf, data) {
  layer(
    data = data,
    mapping = NULL,
    stat = "identity",
    position = "identity",
    inherit.aes = FALSE,                 # <- don't inherit x, y, color, etc.
    geom = ggplot2:::GeomCustomAnn,
    params = list(
      grob = grob,
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    )
  )
}



g_inset1 <- ggplotGrob(g_g1)
g_inset2 <- ggplotGrob(g_g2)
g_inset3 <- ggplotGrob(g_g3)

#add inserts to main plot
g4 <- g_g0 +
  annotation_custom2(
    grob = g_inset1,
    xmin = 0, xmax = 35,
    ymin = 6500, ymax = 9500,
    data = data.frame(SuperGroup = "Eastern")
  )+
  annotation_custom2(
    grob = g_inset2,
    xmin = 0, xmax = 35,
    ymin = 6500, ymax = 9500,
    data = data.frame(SuperGroup = "Fennoscandia")
  )+
  annotation_custom2(
    grob = g_inset3,
    xmin = 0, xmax = 35,
    ymin = 6500, ymax = 9500,
    data = data.frame(SuperGroup = "Southern")
  )

g4

ggsave("Plots/HapNeLD_together_all.pdf", width = 12, height = 7)


# Add inbreeding by age run roh.R
groups_noAsian <- groups[groups$Group %in% order_pop2, ]
groups_noAsian <- groups_noAsian[match(order_pop2, groups_noAsian$Group), ]
roh_time_df_gather_m_noAsian <- roh_time_df_gather_m[roh_time_df_gather_m$SuperGroup!="Asian",]

order_pop2_short <- c("NBAL","SBAL","EEL","CEL",
                      "KAR","EFIN","WFIN","SCAN",
                      "CAR","DINBAL","NWIB","ITA")
roh_time_df_gather_m_noAsian$GroupID <- factor(
  roh_time_df_gather_m_noAsian$GroupID,
  levels = order_pop2_short,
  ordered = TRUE
)


col_vec <- setNames(groups_noAsian$colour, groups_noAsian$GroupID)
xlab_vec <- setNames(
  paste0("<span style='color:", col_vec, "'>", groups_noAsian$GroupID, "</span>"),
  groups_noAsian$GroupID
)



g_age<- ggplot(roh_time_df_gather_m_noAsian,
               aes(GroupID, value*100, color=Group, shape=Status))
g_age <- g_age + 
  geom_point(position=position_jitterdodge(dodge.width = 0.6, jitter.width = 0.6), size=0.7,alpha=0.7)+
  geom_boxplot(outlier.shape = NA, fill=NA )+  facet_wrap(.~key, nrow=1, scales="free_y")+
  theme_classic()+ 
  scale_color_manual(
    values = setNames(groups$colour, groups$Group),  # hex codes mapped to Group
    breaks = groups$Group,                           # ordering of legend
    labels = groups$Group                            # legend labels
  )+  xlab("")+ ylab("% of Genome in ROH")+
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none")+
  scale_x_discrete(labels = xlab_vec)
g_age

g_final <- g4  / g_age + plot_layout(heights = c(1,0.6)) + plot_annotation(tag_levels = 'A')
#ggsave("Plots/HapNeLD_trajectory_rohage.pdf", height = 8, width = 10)
ggsave("Plots/Figure3.pdf", height = 8, width = 10)

# ONLY PURE SAMPLES, >80% OF ONE CLUSTER HAPLONET ------
# Read files recent Ne with HapNe LD pure samples------
ne_files <- unlist(lapply(1:10, function(i) {
  list.files(paste0("Files/HapNe_pure/HapNe_try", i, "/"),
             pattern = "\\.csv$",
             full.names = TRUE)
}))


groups_ne <- sub('\\ld_hapne_estimate.csv$', '', basename(ne_files))
groups_ne <- gsub('P_','', groups_ne)
groups_ne <- gsub('_',' ', groups_ne)

groups_ne <- gsub('NorthWesternIberian','North Western Iberian', groups_ne)
groups_ne <- gsub('Southern Western Asian','South Western Asian', groups_ne)
groups_ne <- gsub('WesternFinland','Western Finland', groups_ne)


try_ne <- sub(".*HapNe_try([0-9]+).*", "\\1", ne_files)

ne_list <- lapply(ne_files,FUN = function(files) {
  read.csv2(files, sep=",")
})

ne_list2 <- vector("list", length(ne_files))
for (i in seq_along(ne_files)) {
  ne_list2[[i]] <- cbind.data.frame(
    Group = groups_ne[i],
    Rep   = as.integer(try_ne[i]),
    ne_list[[i]]
  )
}

# Final dataframe
ne_df <- do.call(rbind, ne_list2)

ne_df$Group <- factor(ne_df$Group, levels = order_pop, ordered = T)
ne_df_m <- merge(ne_df, groups[,c("Group","colour","SuperGroup")], by="Group")


# average all rep
ne_avg <- ne_df_m %>%
  group_by(Group, TIME) %>%
  summarise(
    Q0.5_median = median(Q0.5, na.rm = TRUE),
    Q0.5_low  = quantile(Q0.5, 0.025, na.rm = TRUE),
    Q0.5_high = quantile(Q0.5, 0.975, na.rm = TRUE),
    .groups = "drop"
  )
ne_avg <- merge(ne_avg, groups[,c("Group","colour","SuperGroup")], by="Group")

# PLOT by SUBGROUP without Asians
df_filt <- ne_df_m %>%
  filter(SuperGroup != "Asian", TIME <= 50)

# Only average
lab_df <- ne_avg %>%
  group_by(Group, SuperGroup) %>%
  filter(SuperGroup != "Asian", TIME <= 50)%>%
  ungroup() %>%
  mutate(
    Q0.5 = Q0.5_median   # rename so geom_text_repel works without changes
  )


#Plotting pure -----

g_g0 <- ggplot()
g_g0 <- g_g0 + 
  # Background: all tries in grey
  geom_line(data = df_filt, aes(x = TIME, y = Q0.5, group = interaction(Group, Rep), color=colour),
            #stat="smooth", method = "gam", formula = y ~ s(x, k = 10),
            linewidth = 0.5,
            alpha = 0.5)+
  # Average
  geom_line(data = ne_avg[ne_avg$SuperGroup!="Asian",],aes(x = TIME, y = Q0.5_median, group = Group, color = colour),
            # stat="smooth", method = "gam", formula = y ~ s(x, k = 10), se = FALSE,  
            linewidth = 1.1)  +
  facet_wrap(.~SuperGroup, scales="free_y",nrow=1)+
  xlab(paste0("Generations ago")) + ylab("Ne")+
  scale_color_identity(guide = "legend", labels = groups$Group, breaks = groups$colour) +
  scale_fill_identity(guide = "legend", labels = groups$Group, breaks = groups$colour) +
  guides(color=guide_legend(title="Population"))+
  
  coord_cartesian(clip = "on", ylim = c(0,10000))+
  scale_x_continuous(limits=c(0,50),
                     n.breaks = 10,labels = scales::comma,expand=expansion(mult = c(0, 0.05)))+
  theme_classic()+ theme(legend.position = "none") +
  geom_text_repel(data = lab_df[lab_df$TIME==50,],
                  aes(TIME, Q0.5_median,label = Group, color=colour),
                  size = 3,
                  nudge_x = -1,         # move slightly left from edge
                  nudge_y = -1,         # move slightly left from edge
                  direction = "y",      # only move vertically
                  box.padding = 0.3,
                  segment.color = "grey70",
                  segment.size = 0.3,
                  max.overlaps = Inf,
                  show.legend = FALSE
  )
g_g0

ggsave("Plots/HapNeLd_50g_allsamples_pure.pdf", height = 8, width = 14)
ggsave("Plots/FigS20_pure.pdf", height = 4, width = 9)
