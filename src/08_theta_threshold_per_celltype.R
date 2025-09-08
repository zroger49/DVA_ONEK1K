#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro@gmail.com
# @Description: Determine the optional threshold for each celltype
# @software version: R=4.3.0

## Load libraries 
library(tidyverse)

# Load data
theta1 <- list()
for (rep in 1:100){
  temp_data <- read.csv(paste0("results/azimuth_annotation_results/07_theta_simulation/theta_est_inflation_results_theta1_rep", rep, ".txt"), sep = " ")
  theta1[[rep]] <- temp_data
}

theta1 <- do.call("rbind", theta1)


theta05 <- list()
for (rep in 1:100){
  temp_data <- read.csv(paste0("results/azimuth_annotation_results/07_theta_simulation/theta_est_inflation_results_theta2_rep", rep, ".txt"), sep = " ")
  theta05[[rep]] <- temp_data
}

theta05 <- do.call("rbind", theta05)


theta01 <- list()
for (rep in 1:100){
  temp_data <- read.csv(paste0("results/azimuth_annotation_results/07_theta_simulation/theta_est_inflation_results_theta10_rep", rep, ".txt"), sep = " ")
  theta01[[rep]] <- temp_data
}

theta01 <- do.call("rbind", theta01)


# Get optional threshold values
margin <- 0.05 # 5%
obj <- 0.90 # 90% of the values fall between 5% of the true dispersion

treshold_list <- data.frame("n_cells" = c(), "mean" = c())

for (n_cells in unique(theta1$nr_cell)){
  cat("Testing for ", n_cells, " cells\n")
  for (m in unique(theta1$mu)){
    # # Test for real dispersion = 1
    data1 <- theta1 %>% filter(nr_cell == n_cells & mu == m)
    top_th <- 1 + margin * 1
    bottom_th <- 1 - margin * 1
    percent <- sum(data1$theta_mle > bottom_th & data1$theta_mle < top_th) / nrow(data1)
    print(paste0("For theta", 1," mu ", m, " ", round(percent*100, 2),"% are between ", bottom_th, " and ", top_th))
    if (percent < obj){
      next
    }
    
    # Test for real dispersion = 0.5
    data2 <- theta05  %>% filter(nr_cell == n_cells & mu == m)
    top_th <- 0.5 + margin * 0.5
    bottom_th <- 0.5 - margin * 0.5
    percent <- sum(data2$theta_mle > bottom_th & data2$theta_mle < top_th) / nrow(data2)
    print(paste0("For theta", 0.5," mu ", m, " ", round(percent*100, 2),"% are between ", bottom_th, " and ", top_th))
    if (percent < obj){
      next
    }
    
    
    # Test for real dispersion = 0.01
    data3 <- theta01  %>% filter(nr_cell == n_cells & mu == m)
    top_th <- 0.1 + margin * 0.1
    bottom_th <- 0.1 - margin * 0.1
    percent <- sum(data3$theta_mle > bottom_th & data3$theta_mle < top_th) / nrow(data2)
    print(paste0("For theta", 0.1," mu ", m, " ", round(percent*100, 2),"% are between ", bottom_th, " and ", top_th))
    if (percent < obj){
      next
    }
    
    treshold_list <- rbind(treshold_list, data.frame("n_cells" = n_cells, "mean" = m))
    break
  }
  
}


