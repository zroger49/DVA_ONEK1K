#!/bin/bash
#SBATCH --account=bsc83
#SBATCH --chdir=/gpfs/projects/bsc83/Projects/scRNAseq/rogerioer/immune_system_variability/src/onek1_azimuth/
#SBATCH --cpus-per-task=112
#SBATCH --qos=gp_bscls
#SBATCH --time=06:00:00
#SBATCH --job-name=sctrasnform
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --output=/gpfs/projects/bsc83/Projects/scRNAseq/rogerioer/immune_system_variability/log_dir/onek1_azimuth/outputs/03_generate_mean_variance_matrix.%A_%a.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/scRNAseq/rogerioer/immune_system_variability/log_dir/onek1_azimuth/errors/03_generate_mean_variance_matrix.%A_%a.err
#SBATCH --array=1

module purge; module load impi/2021.10.0 intel/2023.2.0 gsl/2.7.1 mkl/2023.2.0 gcc/13.2.0 greasy/2.2.4.2 R/4.3.0

Rscript 03_generate_mean_variance_matrix.R --cell_type `sed -n ${SLURM_ARRAY_TASK_ID}p 03_generate_mean_variance.params.tab.params.tab | cut -f1`