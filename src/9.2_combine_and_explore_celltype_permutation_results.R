#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro@gmail.com
# @Description: Combine the results across celltypes and interations and analyse them
# @software version: R=4.3.0

library(tidyverse)

########## SETTING DIRECTORIES ############### 
main.dir <- "results/onek1_author_data/"
data.dir <- paste0(main.dir, "15_model_fit/")
out.dir <- paste0(main.dir, "16_model_fit_combined/")

if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

cell_types <- c("B_IN", "B_MEM", "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET", "CD8_S100B", "DC", 
                "Erythrocytes", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma", "Platelets")

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
