#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro@gmail.com
# @Description: Perform enrichement analysis on the Differential Variable genes
# @software version: R=4.3.0



# Load libs ----
library(tidyverse)
source("src/utils/enrichement.R")

# Set dir -----
data.dir <- "results/azimuth_annotation_results/09_explore_DVA_results/"
data.strat.dir <- "results/azimuth_annotation_results/10_explore_DVA_results_sex_stratified/"
results.dir <- "results/azimuth_annotation_results/_enrichement/"

if (!dir.exists(results.dir)){dir.create(results.dir)}

# Load data-----
DVG.filtered.group <- readRDS(paste0(data.dir, "1.DVG.group.filter.rds"))
DVG.filtered.group.strat <- readRDS(paste0(data.strat.dir, "1.DVG.group.filter.rds"))

# Enrichment analysis----

## Enrichment for data not stratified----

enrichement.age <- lapply(names(DVG.filtered.group$age), function(x) {
  up_genes <- DVG.filtered.group$age[[x]] %>% 
    filter(significance_after_group_filter == "sig" & beta  > 0) %>% 
    pull(gene)
  
  down_genes <- DVG.filtered.group$age[[x]] %>% 
    filter(significance_after_group_filter == "sig" & beta  < 0) %>% 
    pull(gene)
  
  
  background <- DVG.filtered.group$age[[x]] %>% pull(gene)
  
  perform_enrichment(up_genes, down_genes, background)
  
}
)

enrichement.sex <- lapply(names(DVG.filtered.group$sex), function(x) {
  up_genes <- DVG.filtered.group$sex[[x]] %>% 
    filter(significance_after_group_filter == "sig" & beta  > 0) %>% 
    pull(gene)
  
  down_genes <- DVG.filtered.group$sex[[x]] %>% 
    filter(significance_after_group_filter == "sig" & beta  < 0) %>% 
    pull(gene)
  
  
  background <- DVG.filtered.group$sex[[x]] %>% pull(gene)
  
  
  perform_enrichment(up_genes, down_genes, background)
}
)

## Enrichment for stratified data


enrichement.M.strat<- lapply(names(DVG.filtered.group.strat$M), function(x) {
  up_genes <- DVG.filtered.group.strat$M[[x]] %>% 
    filter(significance_after_filter == "sig" & beta  > 0) %>% 
    pull(gene)
  
  down_genes <- DVG.filtered.group.strat$M[[x]] %>% 
    filter(significance_after_filter == "sig" & beta  < 0) %>% 
    pull(gene)
  
  
  background <- DVG.filtered.group.strat$M[[x]] %>% pull(gene)
  
  
  perform_enrichment(up_genes, down_genes, background)
  
}
)



enrichement.Fe.strat<- lapply(names(DVG.filtered.group.strat$M), function(x) {
  up_genes <- DVG.filtered.group.strat$M[[x]] %>% 
    filter(significance_after_filter == "sig" & beta  > 0) %>% 
    pull(gene)
  
  down_genes <- DVG.filtered.group.strat$M[[x]] %>% 
    filter(significance_after_filter == "sig" & beta  < 0) %>% 
    pull(gene)
  
  
  background <- DVG.filtered.group.strat$M[[x]] %>% pull(gene)
  
  
  perform_enrichment(up_genes, down_genes, background)
  
}
)

# Save data

names(enrichement.M.strat) <- names(enrichement.Fe.strat) <- names(DVG.filtered.group.strat$M)

saveRDS(enrichement.M.strat, paste0(results.dir, "enrichements_age_stratified_M.rds"))
saveRDS(enrichement.Fe.strat, paste0(results.dir, "enrichements_age_stratified_Fe.rds"))


names(enrichement.age) <- names(enrichement.sex) <- names(DVG.filtered.group$age)

saveRDS(enrichement.age, paste0(results.dir, "enrichements_age.rds"))
saveRDS(enrichement.sex, paste0(results.dir, "enrichements_sex.rds"))

