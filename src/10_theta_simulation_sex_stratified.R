#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro@gmail.com
# @Description: Simulate theta values across different sample sizes and mean expression
# @software version: R=4.3.2


# Load libraries ----
shhh <- suppressPackageStartupMessages

shhh(library(optparse))
shhh(library(glmGamPoi))
shhh(library(fitdistrplus))
shhh(library(dplyr))
shhh(library(moments))

############################## OPTIONS PARSER ###################################### 

option_list = list(
  make_option(c("--i"), action="store", default=NA, type='integer',
              help="simulation replicate"), 
  make_option(c("--theta"), action="store", default=NA, type='numeric',
              help="theta to simulate")
)


opt = parse_args(OptionParser(option_list=option_list))

i <- opt$i
theta <- opt$theta

# Set directories ----
main.dir <- "/gpfs/projects/bsc83/Projects/scRNAseq/rogerioer/immune_system_variability/results/azimuth_annotation_results/"
out.dir <- paste0(main.dir, "10_theta_simulation_stratified/")
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

# Run the simulation----
nr_cell = c(10, 25, 50, 100, 200, 300)
mu_true = c(0.001,0.01,0.05, seq(0.1,0.5,0.1),seq(1,3,0.5), 5, 10, 30)
rpl = c()
rep = i


for(k in 1:length(nr_cell)){
  
  res = data.frame(mu = mu_true, nr_cell = NA, mu_est = NA, var_est = NA, theta_moe = NA, prop_moe_NA = NA, theta_mle = NA,
                   sw_moe = NA, sw_mle = NA, sk_moe = NA, sk_mle = NA)
  
  for(j in 1:length(mu_true)){
    
    # Simulate a matrix with 1000 individuals x k cells per ind on the same gene with Poisson distribution
    m = 1000
    n = nr_cell[k]
    mu = res$mu[j]
    size = 1 / theta 
    data = matrix(NA, nrow = m, ncol = n)
    
    for(i in 1:m){
      # data[i,] = rpois(n, mu)
      data[i,] = rnbinom(n = n, size = size, mu = mu)
    }
    
    fit = glmGamPoi::glm_gp(as.matrix(data), size_factors = FALSE, verbose = F)
    # This is the original estimate from CR-MLE. No further adjustment or shrinkage is added.
    theta_mle = fit$overdispersions
    
    m_est = rowMeans(fit$Mu, na.rm=T)
    v_est = apply(data, 1, function(x)var(x, na.rm=T))
    
    theta_moe = (v_est-m_est)/(m_est^2)
    
    sw_moe = tryCatch(shapiro.test(theta_moe), error = function(e) {a=list();a["p.value"]=NA;return(a)} )
    sw_mle = tryCatch(shapiro.test(theta_mle), error = function(e) {a=list();a["p.value"]=NA;return(a)} )
    
    res[j,] = c(res$mu[j], n, mean(m_est,na.rm=T), mean(v_est,na.rm=T), mean(theta_moe,na.rm=T), prop_moe_NA = sum(is.na(theta_moe)) / m, mean(theta_mle,na.rm=T), sw_moe$p.value, sw_mle$p.value, 
                skewness(theta_moe,na.rm=T), skewness(theta_mle,na.rm=T))
    
    # print(j)
  }
  
  rpl = rbind(rpl, res)
}

rpl$replicate = rep


write.table(rpl,paste0(out.dir, "theta_est_inflation_results_theta", size, "_rep",rep,".txt"), row.names = F, col.names = T, quote = F)