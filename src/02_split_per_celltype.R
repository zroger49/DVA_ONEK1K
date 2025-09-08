#! /usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro
# @Description: Split data per celltype
# @software version: R=4.4

#### OneK1K data set ####
# OneK1K_Angli_SeuratObj.rds sent by Angli

# To fix the Matrix package in Seurat objects
Csparse_validate = 'CsparseMatrix_validate'

############# SETTING LIBRARIES AND TIME ###############

shhh <- suppressPackageStartupMessages
shhh(library(Seurat))
shhh(library(SingleCellExperiment))

print(Sys.time())
start=Sys.time()

#SETTING DIRECTORIES----
data.dir <- "/gpfs/scratch/bsc83/Data/scRNAseq/Yazar2022/sce_data_objects/"
results.dir <- "/gpfs/scratch/bsc83/Data/scRNAseq/Yazar2022/sce_data_objects/"


# Read data ----
seurat_sctransformed <- readRDS(paste0(data.dir, "AllCells_cell_types_SCTransformed.rds"))

print("Data loaded!")
# Subset by celltype and save results

for (ct_name in names(table(seurat_sctransformed@meta.data$cell_type))){

  subset_sctransformed <- subset(x = seurat_sctransformed, subset = cell_type %in% ct_name)
  ct_name <- gsub(" ", "", ct_name)
  saveRDS(subset_sctransformed, paste0(results.dir, ct_name, "_sctransformed.rds"))
}

