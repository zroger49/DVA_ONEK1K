#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro@gmail.com
# @Description: Model the dispersion with linear models. Code adapted Lluis Frontera from Mele-Lab.
# @software version: R=4.3.0

# Load libraries

shhh <- suppressPackageStartupMessages
shhh(library(tidyverse))

# Load function for DV
shhh(source("src/utils/myDE.R"))

############################## OPTIONS PARSER ###################################### 
cell_types <- c("B_IN", "B_MEM", "CD4_ET", "CD4_NC", "CD8_ET", "CD4_SOX4", "CD8_NC", "CD8_S100B", "DC", 
                "Erythrocytes", "Mono_C", "Mono_NC", "NK", "NK_R", "Plasma", "Platelets")



########## SETTING DIRECTORIES ############### 
#main.dir <- "/gpfs/projects/bsc83/Projects/scRNAseq/rogerioer/immune_system_variability/results/onek1_author_data/"
main.dir <- "results/onek1_author_data/"
data.dir <- paste0(main.dir, "14_lognorm/")
out.dir <- paste0(main.dir, "15_model_fit/")

if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}


# Set variables ----
covariates <- "age, sex"

# Model gene expression dispersion model ----

run_models <- function(iter){
  cat("Iteration ", iter, " .... \n")  
  for (ct_name in cell_types){
    print(paste0("Running analysis for ", ct_name))
    
    ## Read dispersion table
    PhxNew <- read.table(paste0(data.dir, ct_name,"_log_nrm_dispersion_", iter, ".tsv"), header = T, check.names = F)
    
    #Set gene as the rowname
    rownames(PhxNew) <- PhxNew$gene
    PhxNew$gene <- NULL
    
    ## Read metadata 
    metadata <- readRDS("data/donor_metadata.rds") %>%
      mutate(sex = as.factor(Gender),
             sample_id = assignment, 
             age = Age) %>% 
      select(sample_id, sex, age)
    
    donors <- intersect(metadata$sample_id, colnames(PhxNew))
    
    metadata.donors <- metadata %>% filter(sample_id %in% donors)
    PhxNew <- PhxNew %>% select(all_of(donors))
    
    cat("Same order between data and metadata", all(colnames(PhxNew) == metadata.donors$sampleid), "\n")
    cat("Processing the analysis with ", length(donors), " donors\n")
    cat("1.1 Differential analysis by donor, genes analyzed:", nrow(PhxNew),"\n")
    
    # Create design matrix
    covariates_vec <- unlist(strsplit(covariates, split = ",")) # Get covariates vector
    form_covs.string <- paste('~', paste0(covariates_vec, collapse = " + "))
    form_covs <- as.formula(form_covs.string) # Get linear model
    
    design_mat <- model.matrix(form_covs, metadata.donors)
    
    rn <- rownames(PhxNew) # Get gene names
    
    ## Start loop by gene
    DVG <- list()
    
    for (phe in c("age", "sex")){
      results <- lapply(rn, function(ii){
        y <- as.numeric(PhxNew[ii, ]) # Get expression vector for gene
        dd <- myDE(y, ii, phenotype = phe, design_mat)
        dd
      })
      
      results <- results[!sapply(results, is.null)] # Filter NULL values
      results <- as.data.frame(do.call(rbind, results)) %>% mutate(celltype = ct_name)
      results$p_value_adj_own <- p.adjust(results$p_value_own,  method = "BH") # Add own adjusted p-value 
      results$p_value_adj <- p.adjust(results$p_value,  method = "BH")# Add adjusted p-value 
      rownames(results) <- NULL
      
      DVG[[phe]] <- results
    }
    
    
    opfn <- paste0(out.dir, iter, ".DVG_results.", ct_name, ".rds")
    saveRDS(DVG, file = opfn)
    print(paste("DVA results for", ct_name, "saved in", opfn))
  }
}

# Apply the function to the list
results <- lapply(10, run_models)

cat ("Finished \n")