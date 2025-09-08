#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro@gmail.com
# @Description: Combine the results across celltypes and interations and analyse them
# @software version: R=4.3.0

library(tidyverse)

########## SETTING DIRECTORIES ############### 
main.dir <- "results/azimuth_annotation_results/"
data.dir <- paste0(main.dir, "17_model_fit/")
out.dir <- paste0(main.dir, "17_model_fit_combined/")

if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}


cell_types <- c(
  "Bintermediate", "Bmemory", "Bnaive", 
  "CD14Mono", "CD16Mono", "CD4CTL", "CD4Naive", 
  "CD4Proliferating", "CD4TCM", "CD4TEM", 
  "CD8Naive", "CD8Proliferating", "CD8TCM", "CD8TEM", 
  "Eryth", "HSPC", "MAIT", 
  "NKProliferating", "NK_CD56bright", "NK", 
  "Plasmablast", "Platelet", "Treg", 
  "cDC2", "dnT", "gdT", "pDC"
)

iter <- 1:10

#### Combine the results #####

res <- list()

# Load and parse data
combine_data <- function(iter){
  for (ct in cell_types){
    data <- readRDS(paste0(data.dir, iter, ".DVG_results.", ct, ".rds"))
    age <- data$age
    sex <- data$sex
    
    
    res[["age"]][[ct]] <- age
    res[["sex"]][[ct]] <- sex
  }
  
  #  Save data ----
  saveRDS(res, paste0(out.dir, iter, ".DVG_results.rds"))
}

lapply(1:10, combine_data)
