#! /usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro
# @Description: Normalize the data based on sequencing depth and pool
# @software version: R=4.4.0

shhh <- suppressPackageStartupMessages
shhh(library(Seurat))


# setup dirs ----
# Main data dirs
data.dir <- "/gpfs/projects/bsc83/Data/scRNAseq/Yazar2022/sce_data_objects/"
# results dir
results.dir <- "/gpfs/scratch/bsc83/Data/scRNAseq/Yazar2022/sce_data_objects/"

#if(!dir.exists(results.dir)){dir.create(results.dir)}


# Load data----
so <- readRDS(paste0(data.dir, "AllCells_cell_type_sceraw.rds"))

# Tranform to Seurat
so <- as.Seurat(so)

# SCTtransform the data
so <- SCTransform(object = so, vars.to.regress = c("date", "percent.mt"), assay = "originalexp")

saveRDS(so, paste0(results.dir, "AllCells_cell_types_SCTransformed.rds"))