#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro@gmail.com
# @Description: Model the dispersion with linear models. Code adapted Lluis Frontera from Mele-Lab
# @software version: R=4.3.0

# Load libraries
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
shhh(library(dplyr))


# Combine results
out.dir <- "results/azimuth_annotation_results/06_model_fit/"

if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
res <- list()

cell_types <- c("CD4CTL",  "CD8TCM", "NKProliferating", "cDC2", "CD4Naive", "CD8TEM", 
                "NK_CD56bright", "dnT", "Bintermediate", "CD4Proliferating", "NK", "gdT", "Bmemory",
                "CD4TCM", "Eryth", "Plasmablast", "pDC", "Bnaive", "CD4TEM", "HSPC", "Platelet", "CD14Mono", 
                "CD8Naive", "Treg", "CD16Mono", "CD8Proliferating", "MAIT")

# Load and parse data
for (ct in cell_types){
  data <- readRDS(paste0(out.dir, "1.DVG_results.", ct, ".rds"))
  age <- data$age
  sex <- data$sex


  res[["age"]][[ct]] <- age
  res[["sex"]][[ct]] <- sex
}

## Save data
saveRDS(res, paste0(out.dir, "1.DVG_results.rds"))

# Combine sex stratified results

# Load and parse data
for (ct in cell_types){
  data <- readRDS(paste0(out.dir, "1.DVG_results.sex_stratified", ct, ".rds"))
  M <- data$M
  Fe <- data$F


  res[["M"]][[ct]] <- M
  res[["F"]][[ct]] <- Fe
}

## Save data
saveRDS(res, paste0(out.dir, "1.DVG_results.sex_stratified.rds"))
