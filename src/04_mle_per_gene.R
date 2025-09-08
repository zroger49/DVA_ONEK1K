#! /usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro
# @Description: Code to compute the dispersion per celltype across genes
# @software version: R=4.4

#### OneK1K data set ####

# To fix the Matrix package in Seurat objects
Csparse_validate = 'CsparseMatrix_validate'

############# SETTING LIBRARIES AND TIME ############### 

shhh <- suppressPackageStartupMessages
shhh(library(optparse))
shhh(library(Seurat))
shhh(library(glmGamPoi))

print(Sys.time())
start=Sys.time()

############################## OPTIONS PARSER ###################################### 

option_list = list(
  make_option(c("--cell_type"), action="store", default=NA, type='character',
              help="Cell types in low resolution (predicted.celltype.l1) or high resolution (cell_type)")
)

opt = parse_args(OptionParser(option_list=option_list))

ct_name <- opt$cell_type

########## SETTING DIRECTORIES ############### 

data.dir <- "/gpfs/scratch/bsc83/Data/scRNAseq/Yazar2022/sce_data_objects/"
main.dir <- "/gpfs/projects/bsc83/Projects/scRNAseq/rogerioer/immune_system_variability/results/azimuth_annotation_results/"

out.dir <- paste0(main.dir, "04_individual_dispersions/")
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}


# Read count matrix ----
print("Reading count matrix")
print(ct_name)
print(data.dir)
ok1k <- readRDS(paste0(data.dir, ct_name, "_sctransformed.rds"))

# Extract a list of individual IDs ----
a <- unlist(table(ok1k@meta.data$assignment))
ind <- as.vector(names(a))

res <- as.data.frame(matrix(NA, nrow=nrow(ok1k[["SCT"]]), ncol=length(ind) + 1)) # + 1 for gene row 
colnames(res) <- c("gene", ind) # Set column names ("gene" and donor)
res$gene = rownames(ok1k[["SCT"]]) # Assign gene names


for (k in 1:length(ind)){
  # Estimate the within-individual dispersion
  
  print(paste("Estimate the within-individual dispersion for donor", k, paste0("(", ind[as.numeric(k)] , ")")))
  
  # For each individual, subset the data set
  ss <- subset(x = ok1k, subset = assignment == as.character(ind[as.numeric(k)])) # Extract data from current individual
  tmp <- GetAssayData(object = ss, assay = "SCT", slot = "counts") # Extract counts from current individual
  tmp2 <- as.data.frame(tmp) # Transform to dataframe
  
  if(ncol(tmp2) < 5){
    res[,k+1] <- NA # Set "NA" if there are less than 5 cells 
  } else {
    fit <-  glmGamPoi::glm_gp(as.matrix(tmp2), size_factors = FALSE, verbose = F) # Parametrization of dispersion
    res[,k+1] <- fit$overdispersions # Add result to dataframe
  }
  
}

# Save result
write.table(res, paste0(out.dir, ct_name, "_within_ind_dispersion.txt"), row.names=T, col.names=T, quote=F)


print("Script ends")
print(Sys.time())
end=Sys.time()
print(difftime(end, start, units = "mins"))


