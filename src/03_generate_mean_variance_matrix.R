#! /usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro
# @Description: Code to compute the mean variance matrix from the scRNA-seq data object
# @software version: R=4.4

#### OneK1K data set ####
# OneK1K_Angli_SeuratObj.rds sent by Angli splited in each celltype

# To fix the Matrix package in Seurat objects
Csparse_validate = 'CsparseMatrix_validate'

############# SETTING LIBRARIES AND TIME ###############

shhh <- suppressPackageStartupMessages
shhh(library(optparse))
shhh(library(Seurat))
shhh(library(SingleCellExperiment))

print(Sys.time())
start=Sys.time()

################# PARSER ARGUMENTS #############

option_list = list(
  make_option(c("--cell_type"), action="store", default=NA, type='character',
              help="Cell types in low resolution (predicted.celltype.l1) or high resolution (cell_type)")
)
opt = parse_args(OptionParser(option_list=option_list))

# Cell type naming for further analysis
ct_name <- opt$cell_type

#SETTING DIRECTORIES----

data.dir <- "/gpfs/scratch/bsc83/Data/scRNAseq/Yazar2022/sce_data_objects/"
main.dir <- "/gpfs/projects/bsc83/Projects/scRNAseq/rogerioer/immune_system_variability/results/azimuth_annotation_results/"

outdir <- paste0(main.dir, "03_mean_variance_estimates/")
if(!dir.exists(outdir)){dir.create(outdir, recursive = T)}


# Read data ----
seurat_sctransformed <- readRDS(paste0(data.dir, ct_name, "_sctransformed.rds"))

print("Data loaded!")

# Extract a list of individual IDs
a <- unlist(table(seurat_sctransformed@meta.data$assignment))
ind <- as.vector(names(a))

## Histogram for cell distribution across individuals
pdf(paste0(outdir, ct_name, "_cells_per_ind.pdf"), width = 6.08, height = 4.59)
hist(a, xlab="Cells per individual", main=paste0("Number of ",ct_name," cells per individual"))
dev.off()


## Check mean-variance relationship within an individual
# To ensure the estimated within-ind variance is accurate, we excluded those individuals with < 5 cells
keep_list <- names(which(a >= 5))
seurat_sctransformed <- subset(x = seurat_sctransformed, subset = assignment %in% keep_list)
print(paste0(length(keep_list)," individuals remain after excluding those with < 5 cells"))


# Use the individual with the most cells as matrix format reference
ok1k_max_ind <- subset(x = seurat_sctransformed, subset = assignment == names(which.max(a)))
ok1k_max_ind <- FindVariableFeatures(object = ok1k_max_ind, selection.method = 'vst', nfeatures = 2000)
hvf.info <- HVFInfo(object = ok1k_max_ind, selection.method = 'vst')

ind <- unique(seurat_sctransformed@meta.data$assignment)
ind <- ind[order(ind)]
gene <- rownames(hvf.info)
print(paste0(length(ind), " individuals and ", length(gene), " genes included in mean and variance matrices!"))

# Save result for individuals and genes analyzed
write.table(as.data.frame(gene),paste0(outdir, ct_name, "_cells_gene_list.txt"), row.names=F, col.names=T, quote=F)

# Create empty matrices
mean_mx <- matrix(NA, nrow=nrow(hvf.info), ncol=length(ind), dimnames = list(gene, ind))
var_mx <- matrix(NA, nrow=nrow(hvf.info), ncol=length(ind), dimnames = list(gene, ind))
mean_mx <- as.data.frame(mean_mx)
var_mx <- as.data.frame(var_mx)

# Checking rownames (genes) and colnames (donors)
print(all(rownames(mean_mx) == rownames(var_mx)) &&
        all(rownames(mean_mx) == rownames(hvf.info)))

print(all(colnames(mean_mx) == colnames(var_mx)) &&
        all(colnames(mean_mx) == ind))

for(i in 1:length(ind)){
  # Subset the data within one individual
  ok1k_ind <- subset(x = seurat_sctransformed, assignment==ind[i])
  ok1k_ind_assay <- ok1k_ind@assays$SCT@data
  
  mean_mx[,i] <- rowMeans(ok1k_ind_assay)
  var_mx[,i] <- rowVars(ok1k_ind_assay)

}

## Save results
write.table(mean_mx, paste0(outdir, ct_name, "_cells_mean_mx.txt"), row.names=T, col.names=T, quote=F)
write.table(var_mx, paste0(outdir, ct_name, "_cells_var_mx.txt"), row.names=T, col.names=T, quote=F)

print("Script ends")
print(Sys.time())
end=Sys.time()
print(difftime(end, start, units = "mins"))
