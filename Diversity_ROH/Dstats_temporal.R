# Ancient Dstats
library(ggplot2)
library(readxl)
library(tidyr)
library(dplyr)
library("patchwork")

setwd("~/Documents/OneDrive - University of Copenhagen/ResQ_MSCA/")

metadata_mw <- read_excel("Files/Wolves_metadata_Juny2026.xlsx", sheet = 1)
samples <- metadata_mw$SampleID

groups <- read_excel("Files/Wolves_metadata_Juny2026.xlsx", sheet = 2)

order_popid <- c("EASIA","VOUR","SWASIA","ZOO","KAR",
                 "EFIN","WFIN", "SCAN", "NBAL",
                 "SBAL","CAR" ,"CEL",
                 "DINBAL","NWIB","ITA")

order_pop <- c("Karelian","Eastern Finland",
               "Western Finland", "Scandinavian", "Northern Baltic",
               "Southern Baltic","Volga-Ural","Central European Lowlands",
               "Carpathian","Dinaric Balkan","North Western Iberian","Italian",
               "Eastern Asian","South Western Asian","Zoo")

order_pop2 <- c( "Karelian","Eastern Finland","Western Finland", "Scandinavian", "Zoo",
                 "Northern Baltic","Southern Baltic","Central European Lowlands",
                 "Volga-Ural" , "Carpathian" ,"Dinaric Balkan","North Western Iberian","Italian")

metadata_mw$Group <- factor(metadata_mw$Group, levels = order_pop, ordered = T)
metadata_mw$GroupID <- factor(metadata_mw$GroupID, levels = order_popid, ordered = T)

groups$Group <- factor(groups$Group, levels = order_pop, ordered = T)
groups$GroupID <- factor(groups$GroupID, levels = order_popid, ordered = T)


metadata_mw2 <- merge(metadata_mw, groups[, c("Group","colour","SuperGroup")], by="Group")

metadata_mw2$SuperGroup <- factor(metadata_mw2$SuperGroup, levels=c("Fennoscandia","Eastern","Southern","Asian","Zoo"), ordered = T)

#metadata ancient samples
ancient_m <- read_excel("Files/filepaths_metadata_edited4.xlsx",sheet = 1)
ancient_m$pop2 <- ancient_m$name


# Read Rdata from Evelyn's Dstats
dstats <- readRDS("Files/Dstat_ancient_SCAN.rds")
dstats$SampleID <- dstats$pop4

dstats_m <- merge(dstats, metadata_mw2[,c("SampleID","Group","colour","SuperGroup")], by="SampleID")
dstats_ma <- merge(dstats_m, ancient_m[,c("pop2","age","age2","Country")], by="pop2")
dstats_ma$ancient <- paste0(dstats_ma$Country,"_", dstats_ma$age)
dstats_ma$age <- as.numeric(dstats_ma$age)

ancient_order <- dstats_ma %>%
  arrange(age) %>%
  distinct(ancient, .keep_all = TRUE) %>%
  pull(ancient)

dstats_ma$ancient <- factor(dstats_ma$ancient, levels = ancient_order, ordered = TRUE)


d_stats_summary <- dstats_ma %>%
  filter(!SuperGroup %in% c("Asian", "Zoo")) %>%
  group_by(SuperGroup,Group,ancient) %>%
  summarise(
    mean_dstats = mean(est, na.rm = TRUE),
    sd_dstats  = sd(est, na.rm = TRUE),   
    n_sign = sum(abs(z) > 3, na.rm = TRUE),
    prop = mean(abs(z) > 3, na.rm = TRUE),
    sign        = ifelse(mean_dstats > 0, "positive", "negative"),
    n_samples = n() ,                      
    .groups = "drop"
  ) %>% arrange(desc(mean_dstats))



d_stats_summary <- d_stats_summary %>%
  mutate(
    signed_prop = ifelse(sign == "negative", -prop, prop)
  )

d_stats_summary2 <- d_stats_summary %>%
  left_join(groups[,c("GroupID","Group","colour")], by = "Group") %>%
  mutate(Group_lab = paste0("<span style='color:", colour, "'>", GroupID, "</span>"))

d_stats_summary2$ancient <- factor(d_stats_summary2$ancient, levels = ancient_order, ordered = TRUE)
d_stats_summary2$Group<- factor(d_stats_summary2$Group, levels=order_pop2, ordered = T)

my_shapes <- c(
  # Group 1: Empty / Hollow Symbols
  0, 1, 2, 5,       # Square, Circle, Triangle up, Diamond
  
  # Group 2: Crosses / Line Symbols
  3, 4, 10, 8,      # Plus, Cross (X), Circle with plus, Star
  
  # Group 3: Full / Filled Symbols
  15, 16, 17, 18    # Filled square, Filled circle, Filled triangle, Filled diamond
)

# plot 
gg <-ggplot(d_stats_summary, aes(x = ancient, y = mean_dstats, shape=Group)) 
gg <- gg+ # geom_errorbar(aes(ymin = mean_dstats - sd_dstats, ymax = mean_dstats + sd_dstats), width = 0.2, color = "grey60") +
  geom_point(aes(color = signed_prop), size = 3) +
  scale_shape_manual(values= my_shapes)+
  scale_color_gradient2(low = "blue", mid = "grey70", high = "red", midpoint = 0, limits = c(-1, 1), name = "Proportion\nsignificant") +
  facet_grid( ~ SuperGroup, scales = "free_x", space = "free_x", drop = TRUE) +
  labs(x = "", y = "Mean D-statistic") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gg

ggsave("Plots/Dstat_test1.pdf", width = 10, height = 7)


gg <-ggplot(d_stats_summary, aes(x = ancient, y = mean_dstats, shape=Group)) 
gg <- gg+ # geom_errorbar(aes(ymin = mean_dstats - sd_dstats, ymax = mean_dstats + sd_dstats), width = 0.2, color = "grey60") +
  geom_point(aes(color = signed_prop), size = 3) +
  scale_shape_manual(values= my_shapes)+
  scale_color_gradient2(low = "blue", mid = "grey70", high = "red", midpoint = 0, limits = c(-1, 1), name = "Proportion\nsignificant") +
    facet_grid(SuperGroup ~., scales = "free_x", space = "free_x", drop = TRUE) +
  labs(x = "", y = "Mean D-statistic") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gg

ggsave("Plots/Dstat_test2.pdf", width =8, height = 10)


#facet plot idea
plots_by_supergroup <- d_stats_summary2 %>%
  split(.$SuperGroup) %>%
  purrr::keep(~ nrow(.x) > 0) %>%
  lapply(function(df) {
    ggplot(df, aes(x = ancient, y = mean_dstats)) +
      geom_errorbar(aes(ymin = mean_dstats - sd_dstats, ymax = mean_dstats + sd_dstats), width = 0.2, color = "grey60") +
      geom_point(aes(color = signed_prop), size = 3) +
      scale_color_gradient2(low = "blue", mid = "grey80", high = "red", midpoint = 0, limits = c(-1, 1), name = "Proportion\nsignificant") +
      #facet_grid(Group ~ SuperGroup, drop = TRUE)+
      ylim(c(-0.02,0.05))+
      facet_grid(Group_lab ~SuperGroup, scales = "free",  drop = TRUE, switch = "y") +
      labs(x = "", y = "Mean D-statistic") +
      theme_classic() +
      theme(  strip.background.y = element_blank(), legend.position = "none",
        strip.text.y.left = element_markdown(angle = 0),
        strip.placement = "outside",
        axis.text.x = element_text(angle = 45, hjust = 1)
      )  })

plots_by_supergroup[["Fennoscandia"]]
plots_by_supergroup[["Eastern"]]
plots_by_supergroup[["Southern"]]


library(patchwork)

n <- length(plots_by_supergroup)

plots_by_supergroup <- lapply(seq_along(plots_by_supergroup), function(i) {
  if (i == n) {
    plots_by_supergroup[[i]] + theme(legend.position = "right")
  } else {
    plots_by_supergroup[[i]] + theme(legend.position = "none")
  }
})

patchwork::wrap_plots(plots_by_supergroup, nrow = 1)



#facet plot idea but with shape
plots_by_supergroup <- d_stats_summary2 %>%
  split(.$SuperGroup) %>%
  purrr::keep(~ nrow(.x) > 0) %>%
  lapply(function(df) {
    ggplot(df, aes(x = ancient, y = mean_dstats, shape=Group)) +
      #geom_errorbar(aes(ymin = mean_dstats - sd_dstats, ymax = mean_dstats + sd_dstats), width = 0.2, color = "grey60") +
      geom_point(aes(color = signed_prop), size = 3) +
      scale_color_gradient2(low = "blue", mid = "grey80", high = "red", midpoint = 0, limits = c(-1, 1), name = "Proportion\nsignificant") +
      ylim(c(-0.02,0.05))+
      facet_grid(.~SuperGroup, scales = "free",  drop = TRUE, switch = "y") +
      labs(x = "", y = "Mean D-statistic") + guides(color="none")+
      theme_classic() +
      theme(  legend.position = c(0.2,0.9),
              legend.title=element_blank(),
              strip.placement = "outside",
              axis.text.x = element_text(angle = 45, hjust = 1)
      )  })

plots_by_supergroup[["Fennoscandia"]]
plots_by_supergroup[["Eastern"]]
plots_by_supergroup[["Southern"]]


library(patchwork)

n <- length(plots_by_supergroup)

plots_by_supergroup <- lapply(seq_along(plots_by_supergroup), function(i) {
  if (i == n) {
    plots_by_supergroup[[i]] + theme(legend.position = "right")
  } else {
    plots_by_supergroup[[i]] + theme(legend.position = "none")
  }
})

patchwork::wrap_plots(plots_by_supergroup, nrow = 1)




