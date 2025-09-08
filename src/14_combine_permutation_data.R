#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro@gmail.com
# @Description: Combine the data from the permutation analaysis
# @software version: R=4.3.0

# Load libraries

shhh <- suppressPackageStartupMessages
shhh(library(tidyverse))
shhh(library(optparse))
shhh(library(ggpubr))

########## SETTING DIRECTORIES ############### 
main.dir <- "results/azimuth_annotation_results/"
data.dir <- paste0(main.dir, "13_model_fit_permutations/")
out.dir <- paste0(main.dir, "14_model_fit_permutations_analysis/")

if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

# Celltypes 
#cell_types <- c("Bintermediate", "Bmemory", "Bnaive", "CD4CTL", "CD4Naive", "CD4TCM", "CD4TEM", "CD8Naive",
#                  "CD8TCM", "CD8TEM", "CD14Mono", "CD16Mono","MAIT", "NK_CD56bright", "NK", "Treg")

cell_types <- c("CD4Naive", "CD4TEM", "CD8TEM", "NK")


# Load permutation data summary data

summary <- list()
#analysis <- list()
for (phe in c("sex", "age")){
  summary[[phe]] <- list()
  #analysis[[phe]] <- list()
  for (ct_name in cell_types){
    #analysis[[phe]][[ct_name]] <- list()
    ct_phe <- data.frame()
    for (i in 1:10){
      data <- readRDS(paste0(data.dir, "DVG.permutation.iter", i, ".", phe, ".", ct_name, ".summary.rds"))
      data$ct <- ct_name
      ct_phe <- rbind(ct_phe, data)
    }
    summary[[phe]][[ct_name]] <- ct_phe
  }
}


# Combine all the datasets into a single table
combined_data.age <- do.call("rbind", summary$age)
combined_data.sex <- do.call("rbind", summary$sex)


# Load and combine with the analysis data
#DVG_after_filter <- readRDS("results/onek1_author_data/main_results_exploration/DVG.group.filter.rds")
#main.dir
DVG_after_filter <- readRDS(paste0(main.dir, "/09_explore_DVA_results/1.DVG.group.filter.rds")) 

real_data.sex <- do.call("c", lapply(DVG_after_filter$sex, function(x) sum(x$significance_after_group_filter == "sig")))
real_data.age <- do.call("c", lapply(DVG_after_filter$age, function(x) sum(x$significance_after_group_filter == "sig")))

real_data.sex <- real_data.sex[names(real_data.sex) %in% unique(combined_data.age$ct)]
real_data.age <- real_data.age[names(real_data.age) %in% unique(combined_data.age$ct)]


min_pvalues_after_filtering.sex <-  do.call("c", lapply(DVG_after_filter$sex, function(x) min(x$p_value, na.rm = TRUE)))
min_pvalues_after_filtering.age <-  do.call("c", lapply(DVG_after_filter$age, function(x) min(x$p_value, na.rm = TRUE)))

min_pvalues_after_filtering.sex <- min_pvalues_after_filtering.sex[names(min_pvalues_after_filtering.sex) %in% unique(combined_data.age$ct)]
min_pvalues_after_filtering.age <- min_pvalues_after_filtering.age[names(min_pvalues_after_filtering.age) %in% unique(combined_data.age$ct)]



real_data.sex <- data.frame("iter" = "real", "phe" = "sex", min_value = NA, 
                            min_pval_after_filter = min_pvalues_after_filtering.sex, number_of_DVG_before_filter  = NA, 
                            number_of_DVG_after_filter = real_data.sex, permutation_metric = NA, ct = names(real_data.sex))


real_data.age <- data.frame("iter" = "real", "phe" = "age", min_value = NA, 
                            min_pval_after_filter = min_pvalues_after_filtering.age, number_of_DVG_before_filter  = NA, 
                            number_of_DVG_after_filter = real_data.age, permutation_metric = NA, ct = names(real_data.age))

combined_data.age.real <- rbind(combined_data.age, real_data.age)
combined_data.sex.real <- rbind(combined_data.sex, real_data.sex)


combined_data.age.real$status <- ifelse(combined_data.age.real$iter == "real", "real", "simulation")
combined_data.sex.real$status <- ifelse(combined_data.sex.real$iter == "real", "real", "simulation")



a <- ggplot(data = combined_data.age.real, aes(x = ct, y = number_of_DVG_after_filter, colour = status)) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.7) +
  labs(
    title = "Scatterplot (Age)",
    color = " "
  ) +
  scale_colour_manual(values = c("simulation" = "blue", "real" = "darkred")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    axis.text = element_text(angle = 90, vjust = 1)
  ) +xlab("") + ylab("")


b <- ggplot(data = combined_data.sex.real, aes(x = ct, y = number_of_DVG_after_filter, colour = status)) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.7) +
  labs(
    title = "Scatterplot (sex)",
    color = " "
  ) +
  scale_colour_manual(values = c("simulation" = "blue", "real" = "darkred")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    axis.text = element_text(angle = 90, vjust = 1)
  ) +xlab("") + ylab("")

pdf(paste0(out.dir, "permutation_analysis.pdf"), w  = 8, h = 6)
ggarrange(a,b, ncol = 1, nrow = 2, common.legend = TRUE)
dev.off()
