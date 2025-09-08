#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro@gmail.com
# @Description: Prepare the dispersion estimates for modeling. Prepare one dispersion matrices per iteration and per cell type
# @software version: R=4.3.2

shhh <- suppressPackageStartupMessages
shhh(library(optparse))
shhh(library(ggplot2))
shhh(library(gridExtra))
shhh(library(data.table))
shhh(library(R.utils))

############################## OPTIONS PARSER ###################################### 

option_list = list(
  make_option(c("--cell_type"), action="store", default=NA, type='character',
              help="Cell types in low resolution (predicted.celltype.l1) or high resolution (cell_type)"),
  make_option(c("--iteration"), action="store", default=NA, type='character',
              help="Number of iteration")
)

opt = parse_args(OptionParser(option_list=option_list))

ct_name <- opt$cell_type
iter <- as.numeric(opt$iteration)

########## SETTING DIRECTORIES ############### 

main.dir <- "/gpfs/projects/bsc83/Projects/scRNAseq/rogerioer/immune_system_variability/results/azimuth_annotation_results/"

out.dir <- paste0(main.dir, "16_lognorm/")
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

img_dir <- paste0(out.dir, "images/")
if(!dir.exists(img_dir)){
  dir.create(img_dir)
}

##################

#Read pseudo-bulk matrix (mean) 
#-> This is required to remove genes no expressed 
#in over 10% of the individuals.
mx <- fread(paste0(main.dir, "03_mean_variance_estimates/", ct_name, "_cells_mean_mx.txt"), header = F)
mx <- as.data.frame(mx)
mx <- mx[, -1] #Remove gene annotation

# Remove individuals with all NAs
index_filter <- FALSE
if(sum(is.na(mx[1, ]))>0){
  index_filter = TRUE
  index = which(is.na(mx[1, ]))
  mx = mx[ ,-index]
}
pi0 = rowSums(mx==0)/ncol(mx)

# Read pseudo-bulk matrix (dispersion)
mxd <- read.csv(paste0(main.dir, "15_individual_dispersions_cell_permutations/", ct_name, "_within_ind_dispersion_", iter, ".txt"), header = T, sep = " ")
mxd = as.data.frame(mxd)
#Remove the X from the colnames
colnames(mxd) <- gsub(pattern = "X", replacement = "", x = colnames(mxd))
# Remove genes names for now
gene_names <- mxd[,1]
mxd = mxd[,-1]

cat(ct_name, " has ", nrow(mxd), " genes across ", ncol(mx), "individuals... \n")

if (index_filter){
  mxd = mxd[, -index]
}

# Remove genes with pi0 > 0.9
mxd = mxd[which(pi0 <= 0.9), ]
gene_names <- gene_names[which(pi0 <= 0.9)]


# log(x+1) and standardization
mxd = apply(mxd, 2, function(x)log(x+1))
mxd = as.data.frame(t(scale(t(mxd))))

cat("After filtering...", nrow(mxd), " genes across ", ncol(mx), "individuals... \n")

#Plot data distribution in 5 example genes
set.seed(123)  # Optional: for reproducibility
random_numbers <- sample(1:nrow(mxd), 6, replace = FALSE)

plot_list <- list()  
for (i in random_numbers) {
  # Extract the ith row of mxd and convert it to a numeric vector
  data_i <- as.numeric(mxd[i,])
  
  # Create a histogram using ggplot
  p <- ggplot(data = data.frame(value = data_i), aes(x = value)) +
    geom_histogram(binwidth = 0.01, fill = 'blue', color = 'black') +
    ggtitle(gene_names[i]) +
    theme_minimal()
  
  # Append the plot to the plot_list
  plot_list[[length(plot_list) + 1]] <- p
}

pdf(paste0(img_dir, ct_name, iter, "_random_genes_distribution.pdf"), h = 8,  w = 8)
do.call(grid.arrange, plot_list)
dev.off()

# Add back gene information
mxd$genes <- gene_names

# Save data
write.table(mxd, paste0(out.dir, ct_name, "_log_nrm_dispersion_", iter, ".tsv"), row.names = F, col.names = T, quote = F, sep = "\t")
