#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro@gmail.com
# @Description: Model the dispersion with linear models. Code adapted Lluis Frontera from Mele-Lab
# @software version: R=4.3.0

# Load libraries

shhh <- suppressPackageStartupMessages
shhh(library(optparse))
shhh(library(dplyr))


# Load function for DV
shhh(source("../utils/myDE.R")) # Compatatible with the .sh file 


############################## OPTIONS PARSER ###################################### 

option_list = list(
  make_option(c("--cell_type"), action="store", default=NA, type='character',
              help="Cell types in low resolution (predicted.celltype.l1) or high resolution (cell_type)")
)

opt = parse_args(OptionParser(option_list=option_list))

ct_name <- opt$cell_type

########## SETTING DIRECTORIES ############### 
main.dir <- "/gpfs/projects/bsc83/Projects/scRNAseq/rogerioer/immune_system_variability/results/azimuth_annotation_results/"
#main.dir <- "results/azimuth_annotation_results//"
data.dir <- paste0(main.dir, "05_lognorm/")
out.dir <- paste0(main.dir, "06_model_fit/")

if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

# Age and Sex analysis 

# Set variables ----
covariates <- "age, sex"

# Model gene expression dispersion model ----
print(paste0("Running analysis for ", ct_name))

## Read dispersion table
PhxNew <- read.table(paste0(data.dir, ct_name,"_log_nrm_dispersion.tsv"), header = T, check.names = F)

#Set gene as the rowname
rownames(PhxNew) <- PhxNew$gene
PhxNew$gene <- NULL

## Read metadata 
metadata_file <- "/gpfs/projects/bsc83/Projects/scRNAseq/rogerioer/immune_system_variability/data/donor_metadata.rds"
metadata <- readRDS(metadata_file) %>%
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


opfn <- paste0(out.dir, "1.DVG_results.", ct_name, ".rds")
saveRDS(DVG, file = opfn)
print(paste("DVA results for", ct_name, "saved in", opfn))

print("Running stratified models")

# Age, stratified by age

# Set variables ----
covariates <- "age, sex"

# Model gene expression dispersion model ----

print(paste0("Running analysis for ", ct_name))

## Read dispersion table
PhxNew <- read.table(paste0(data.dir, ct_name,"_log_nrm_dispersion.tsv"), header = T, check.names = F)

#Set gene as the rowname
rownames(PhxNew) <- PhxNew$gene
PhxNew$gene <- NULL

## Read metadata 
metadata_file <- "/gpfs/projects/bsc83/Projects/scRNAseq/rogerioer/immune_system_variability/data/donor_metadata.rds"
metadata <- readRDS(metadata_file) %>%
  mutate(sex = as.factor(Gender),
         sample_id = assignment, 
         age = Age) %>% 
  select(sample_id, sex, age)

donors_all <- intersect(metadata$sample_id, colnames(PhxNew))


## Start loop by gene
DVG <- list()

for (s in c("M", "F")){
  
  metadata.donors <- metadata %>% filter(sample_id %in% donors_all) %>% filter(sex == s)
  donors <- donors_all[donors_all %in% metadata.donors$sample_id]
  PhxNew_sex <- PhxNew %>% select(all_of(donors))
  
  cat("Same order between data and metadata", all(colnames(PhxNew_sex) == metadata.donors$sampleid), "\n")
  cat("Processing the analysis with ", length(donors), " donors\n")
  cat("1.1 Differential analysis by donor, genes analyzed:", nrow(PhxNew_sex),"\n")
  
  # Create design matrix
  covariates_vec <- unlist(strsplit(covariates, split = ",")) # Get covariates vector
  form_covs.string <- paste('~', paste0(covariates_vec, collapse = " + "))
  form_covs <- as.formula(form_covs.string) # Get linear model
  
  design_mat <- model.matrix(form_covs, metadata.donors)
  
  rn <- rownames(PhxNew_sex) # Get gene names
  
  results <- lapply(rn, function(ii){
    y <- as.numeric(PhxNew_sex[ii, ]) # Get expression vector for gene
    dd <- myDE(y, ii, phenotype = "age", design_mat)
    dd
  })
  
  results <- results[!sapply(results, is.null)] # Filter NULL values
  results <- as.data.frame(do.call(rbind, results)) %>% mutate(celltype = ct_name)
  results$p_value_adj_own <- p.adjust(results$p_value_own,  method = "BH") # Add own adjusted p-value 
  results$p_value_adj <- p.adjust(results$p_value,  method = "BH")# Add adjusted p-value 
  rownames(results) <- NULL
  
  DVG[[s]] <- results
}


opfn <- paste0(out.dir, "1.DVG_results.sex_stratified", ct_name, ".rds")
saveRDS(DVG, file = opfn)
print(paste("DVA results for", ct_name, "saved in", opfn))

cat ("Finished \n")
