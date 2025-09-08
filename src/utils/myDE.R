library(lme4)
library(lmerTest)


# Set functions to use ----
myDE <- function(y, gene, phenotype, design_matrix, threshold = 3, random_effects = FALSE) {
  # DEA for a single contrast (batch together)
  # y: Vector of dispersion values
  # gene: Gene identifier
  # design_matrix: Design matrix for covariates
  # threshold: A threshold value for the number of samples (default is 3) 
  if (random_effects == TRUE){
    lm0 <- try(lmer(design_matrix, data = y), silent = TRUE)
  }else{
    lm0 <- try(lm(y ~ 0 + design_matrix), silent = TRUE) # Fit linear model
  }
  summary_fit <- summary(lm0)
  
  dd_list <- list()
  # Extract summary statistics from model
  if (class(y) != "try-error" && class(lm0) != "try-error") {
    if (random_effects == TRUE){
      b <- summary_fit[["coefficients"]]
      coef_interest <- phe
      b <- b[coef_interest, "Estimate"]
    }else{
      b <- coef(lm0) # Get coefficients
      if(phenotype == "age"){
        coef_interest <- "design_matrixage"
      } else {
        coef_interest <- "design_matrixsexF"
      }
      b <- b[coef_interest]
    }
    

    vb <- diag(vcov(lm0)) # Get diagonal of covariance matrix
    vb <- vb[coef_interest]
    std_error <- summary_fit$coefficients[, "Std. Error"]  # Get std error
    std_error <- std_error[coef_interest]
    z_score <- summary_fit$coefficients[, "t value"]  # Get t value
    z_score <- z_score[coef_interest] 
    p_value <- summary_fit$coefficients[, "Pr(>|t|)"]  # Get p-value
    p_value <- p_value[coef_interest] # Get p-value of coefficient of interest
    
    # Extract own summary statistics
    bhat <- b
    sdhat <- sqrt(vb) 
    z <- bhat / sdhat
    p <- 2 * pnorm(-abs(z))
    dd <- data.frame(gene = gene, beta = bhat, stderr_own = sdhat, z_score_own = z, p_value_own = p, 
                     stderr = std_error, z_score = z_score, p_value = p_value) # Merge all into dataframe 
    
  } else {
    dd <- data.frame(gene = gene, beta = NA, stderr_own = NA, z_score_own = NA, p_value_own = NA, 
                     stderr = NA, z_score = NA, p_value = NA)
  }
  dd
}



myDE_stimulation <- function(y, gene, design_matrix, contrast, threshold = 3, random_effects = FALSE) {
  # DEA for a single contrast (batch together)
  # y: Vector of dispersion values
  # gene: Gene identifier
  # design_matrix: Design matrix for covariates
  # threshold: A threshold value for the number of samples (default is 3) 
  # Contrast: The constrat to test
  
  
  # The difference from this function to the previous is that in this function is for the "stimulation" case. 
  # I divided because in this case I need to specify contrasts
  if (random_effects == TRUE){
    lm0 <- try(lmer(design_matrix, data = y), silent = TRUE)
  }else{
    lm0 <- try(lm(y ~ 0 + design_matrix), silent = TRUE) # Fit linear model
  }
  summary_fit <- summary(lm0)
  
  dd_list <- list()
  # Extract summary statistics from model
  if (class(lm0) != "try-error"){
  #if (class(y) != "try-error" && class(lm0) != "try-error") {
    if (random_effects == TRUE){
      b <- summary_fit[["coefficients"]]
      coef_interest <- contrast
      b <- b[coef_interest, "Estimate"]
    }else{
      b <- coef(lm0) # Get coefficients
      coef_interest <- paste0("design_matrix", contrast)
      b <- b[coef_interest]
    }
    
    
    vb <- diag(vcov(lm0)) # Get diagonal of covariance matrix
    vb <- vb[coef_interest]
    std_error <- summary_fit$coefficients[, "Std. Error"]  # Get std error
    std_error <- std_error[coef_interest]
    z_score <- summary_fit$coefficients[, "t value"]  # Get t value
    z_score <- z_score[coef_interest] 
    p_value <- summary_fit$coefficients[, "Pr(>|t|)"]  # Get p-value
    p_value <- p_value[coef_interest] # Get p-value of coefficient of interest
    
    # Extract own summary statistics
    bhat <- b
    sdhat <- sqrt(vb) 
    z <- bhat / sdhat
    p <- 2 * pnorm(-abs(z))
    # Extract error message
    if (random_effects == T){
      msg <- lm0@optinfo$conv$lme4$messages # #"boundary (singular) fit: see help('isSingular')"
      if (is.null(msg)){
        msg <- NA 
      }
      
      if (length(msg) >= 2){
        msg <- paste0(msg, collapse = ";")
      }  
    }else{
      msg <- NA
    }
    
    
    dd <- data.frame(gene = gene, beta = bhat, stderr_own = sdhat, z_score_own = z, p_value_own = p, 
                     stderr = std_error, z_score = z_score, p_value = p_value, msg = msg) # Merge all into dataframe 
    
  } else {
    dd <- NA
  }
  dd
}