#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro@gmail.com
# @Description: Model the dispersion with linear models. Permute the labels and check if there is still significant results
# @software version: R=4.3.0

# Load libraries

shhh <- suppressPackageStartupMessages
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(optparse))

# Load function for DV
shhh(source("../utils/myDE.R"))

############################## OPTIONS PARSER ######################################

option_list = list(
  make_option(c("--cell_type"), action="store", default=NA, type='character',
              help="Cell types in low resolution (predicted.celltype.l1) or high resolution (cell_type)"),
  make_option(c("--phenotype"), action="store", default=NA, type='character',
              help="Phenotype to analyze, either 'sex' or 'age'"),
  make_option(c("--iteration"), action="store", default=NA, type='character',
              help="Number of iterations (e.g., 1 to 10 or more)")
)

# Parsing options
opt <- parse_args(OptionParser(option_list = option_list))

# Cell type naming for further analysis
ct_name <- opt$cell_type
iter <- as.numeric(opt$iteration)
phe <- opt$phenotype


########## SETTING DIRECTORIES ###############
main.dir <- "/gpfs/projects/bsc83/Projects/scRNAseq/rogerioer/immune_system_variability/results/azimuth_annotation_results/"
#main.dir <- "results/onek1_author_data/"
data.dir <- paste0(main.dir, "05_lognorm/")
out.dir <- paste0(main.dir, "13_model_fit_permutations/")

if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

# Determine which genes pass the threshold in each condition
cts_to_study <- c("Bintermediate", "Bmemory", "Bnaive", "CD4CTL", "CD4Naive", "CD4TCM", "CD4TEM", "CD8Naive",
                  "CD8TCM", "CD8TEM", "CD14Mono", "CD16Mono","MAIT", "NK_CD56bright", "NK", "Treg")
filter_mean_ct <- c(10, 10,5,30,1.5,2,10,5,30,2,10,30,30,30,2,10)
names(filter_mean_ct) <- cts_to_study


ct_treshold <- filter_mean_ct[ct_name]

## Load mean expression
mean_exps <- read.table(paste0(main.dir, "03_mean_variance_estimates/", ct_name, "_cells_mean_mx.txt"))

## Filter
### First general filter (per ct_name)
mean_per_gene <- rowMeans(mean_exps, na.rm = TRUE)
keep_genes_first_filter <- names(mean_per_gene[mean_per_gene >= ct_treshold])


### Second filter (group filter)
metadata <- readRDS("../../data/donor_metadata.rds") %>%
  mutate(sex = as.factor(Gender),
         sample_id = assignment,
         age = Age) %>%
  select(sample_id, sex, age)

metadata_to_add <- metadata %>%
  mutate(age_group = ifelse(age < 40, "Y", "O")) %>%
  select(sample_id, sex, age_group)

mean_per_gene <- rowMeans(mean_exps, na.rm = TRUE)

long_df <- mean_exps %>%
  mutate(gene = row.names(.)) %>%  # Add gene row numbers
  pivot_longer(cols = 1:ncol(mean_exps), names_to = "sample_id", values_to = "mean") %>%
  mutate(sample_id = gsub("X", "", sample_id)) %>%
  merge(metadata_to_add, by = "sample_id")

if (phe == "age"){
  group_mean_age <- long_df %>%
    mutate(age_group = as.factor(age_group)) %>%
    group_by(age_group, gene) %>%
    summarise(group_mean = mean(mean))
  
  group_mean_age_filtered <- group_mean_age %>%
    filter(group_mean > ct_treshold)
  
  keep_genes_group_filter.age <- names(table(group_mean_age_filtered$gene)[table(group_mean_age_filtered$gene) == 2])
}else if (phe == "sex"){
  group_mean_sex <- long_df %>%
    mutate(sex = as.factor(sex)) %>%
    group_by(sex, gene) %>%
    summarise(group_mean = mean(mean))
  
  
  group_mean_sex_filtered <- group_mean_sex %>%
    filter(group_mean > ct_treshold)
  
  keep_genes_group_filter.sex <- names(table(group_mean_sex_filtered$gene)[table(group_mean_sex_filtered$gene) == 2])
}

# Set variables ----
covariates <- "age, sex"

analysis <- list()

# Model gene expression dispersion model ----
#for (ct_name in cell_types){
print(paste0("Running analysis for ", ct_name))

## Read dispersion table
PhxNew <- read.table(paste0(data.dir, ct_name,"_log_nrm_dispersion.tsv"), header = T, check.names = F)

#Set gene as the rowname
rownames(PhxNew) <- PhxNew$gene
PhxNew$gene <- NULL

## Read metadata
metadata <- readRDS("../../data/donor_metadata.rds") %>%
  mutate(sex = as.factor(Gender),
         sample_id = assignment,
         age = Age) %>%
  select(sample_id, sex, age)


# I need to do this step before permuting, since some samples might be removed
donors <- intersect(metadata$sample_id, colnames(PhxNew))
metadata.donors <- metadata %>% filter(sample_id %in% donors)

# Permute the labels
metadata.donors$sample_id_permuted <- sample(metadata.donors$sample_id)
donors_new_order <- intersect(metadata.donors$sample_id_permuted, colnames(PhxNew))

PhxNew <- PhxNew %>% select(all_of(donors_new_order))


cat("Same order between data and metadata", all(colnames(PhxNew) == metadata.donors$sample_id_permuted), "\n")
cat("Processing the analysis with ", length(donors), " donors\n")
cat("1.1 Differential analysis by donor, genes analyzed:", nrow(PhxNew),"\n")

# Create design matrix
covariates_vec <- unlist(strsplit(covariates, split = ",")) # Get covariates vector
form_covs.string <- paste('~', paste0(covariates_vec, collapse = " + "))
form_covs <- as.formula(form_covs.string) # Get linear model

design_mat <- model.matrix(form_covs, metadata.donors)

rn <- rownames(PhxNew) # Get gene names

# Get the "distance" metric. This tells me how different the permutation is from the real thing
permuted_labels <- real_labels <- metadata.donors[,phe]
names(permuted_labels) <-  metadata.donors$sample_id_permuted
names(real_labels) <-  metadata.donors$sample_id
real_labels <- real_labels[names(permuted_labels)]

if (phe == "sex"){
  metric <- round(sum(real_labels == permuted_labels) / length(permuted_labels), 3) * 100
}else if (phe == "age"){
  metric <- mean(abs(real_labels - permuted_labels))
}

# Analysis (This is run per gene)
results <- lapply(rn, function(ii){
  y <- as.numeric(PhxNew[ii, ]) # Get expression vector for gene
  dd <- myDE(y, ii, phenotype = phe, design_mat)
  dd
})

results <- results[!sapply(results, is.null)] # Filter NA values
results <- as.data.frame(do.call(rbind, results)) %>% mutate(celltype = ct_name)
results$p_value_adj_own <- p.adjust(results$p_value_own,  method = "BH") # Add own adjusted p-value
results$p_value_adj <- p.adjust(results$p_value,  method = "BH")# Add adjusted p-value
rownames(results) <- NULL
DVG <- results


# Filter DVG based on min expression
if (phe == "age"){
  genes_to_keep_ct <- intersect(keep_genes_first_filter, keep_genes_group_filter.age)
}else if (phe == "sex"){
  genes_to_keep_ct <- intersect(keep_genes_first_filter, keep_genes_group_filter.sex)
}

DVG_filter <- DVG %>%
  filter(gene %in% genes_to_keep_ct) %>%
  mutate(p_value_adj = p.adjust(p_value, method = "BH"))

analysis[["raw"]] <- DVG
analysis[["filtered"]] <- DVG_filter

number_of_DVG_before_filtering <-  sum(DVG$p_value_adj < 0.05)
number_of_DVG_after <-  sum(DVG_filter$p_value_adj < 0.05)
min_pval <- min(DVG$p_value, na.rm = T)
min_pval_after_filter <- min(DVG_filter$p_value, na.rm = T)

summary <- data.frame("iter" = iter, "phe" = phe,
                      "min_value" = min_pval,
                      "min_pval_after_filter" = min_pval_after_filter,
                      "number_of_DVG_before_filter" = number_of_DVG_before_filtering,
                      "number_of_DVG_after_filter" = number_of_DVG_after,
                      "permutation_metric" = metric
)


saveRDS(analysis, paste0(out.dir,"DVG.permutation.iter", iter, ".", phe, ".", ct_name,".resuls.rds"))
saveRDS(summary, paste0(out.dir,"DVG.permutation.iter", iter, ".", phe, ".", ct_name,".summary.rds"))

cat ("Finished \n")