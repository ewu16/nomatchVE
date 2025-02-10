#' Compute simultaneous confidence intervals for timepoints object
#'
#' @param object Object created by call to [timepoints()]
#' @param alpha Significance level for confidence interval
#'
#' @return A dataframe with lower and upper simultaneous confidence intervals for 
#' each timepoint 
#' @export
#'
simultaneous_ci <- function(object, alpha){
  #TODO: check bootstrap samples exists 
  #TODO: add args for specifying time range
  
  #create bootstrap matrix for each term: 
  # rows are bootstrap samples, col are timepoints
  boot_mat <- matrix(NA, nrow = object$n_boot, ncol = length(object$times),
                     dimnames = list(NULL, object$times))
  risk_0_mat <- risk_1_mat <- ve_mat <- boot_mat
  for(i in seq_along(object$times)){
    risk_0_mat[,i] <- object$boot_samples[[i]][,1]
    risk_1_mat[,i] <- object$boot_samples[[i]][,2]
    ve_mat[,i]     <- object$boot_samples[[i]][,3]
  }
  
  #compute z_star for each term
  boot_mat_list <- list(risk_0_mat, risk_1_mat, ve_mat)
  z_star <- sapply(boot_mat_list, \(mat) get_z_star(mat, alpha = alpha))
  
  #compute ci for each timepoint
  simul_ci_list <- lapply(seq_along(object$times), \(i){
    ci <- compute_wald_ci(pt_est = object$estimates[[i]][,"estimate"], 
                    boot_samples = object$boot_samples[[i]],
                    alpha = NULL, z_star = z_star)
    ci_df <- obsve$estimates_to_df(ci)
    ci_df$t0 <- object$times[i]
    ci_df
  })
  
  simul_ci <- do.call("rbind", simul_ci_list)
  names(simul_ci)[names(simul_ci) == "wald_lower"] <- "simul_lower"
  names(simul_ci)[names(simul_ci) == "wald_upper"] <- "simul_upper"
  simul_ci[, c("t0", "term", "simul_lower", "simul_upper")]
}

## get critical value for simultaneous CI
get_z_star <- function(mat, alpha = .05){
  covariance <- tryCatch(stats::cov(mat, use = "complete.obs"),
                         error = function(e){print(e); NaN})
  if(is.null(covariance)){
    return(NaN)
  }
  sds <- sqrt(diag(covariance))
  estimate_draws <- MASS::mvrnorm(n = 10000, mu = rep(0, ncol(mat)), Sigma = covariance)
  max_z <- apply(estimate_draws, 1, \(x) max(x/sds, na.rm = TRUE))
  z_star <- stats::quantile(max_z, 1 - alpha)
  unname(z_star)
}


