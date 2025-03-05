#' Compute simultaneous confidence intervals for timepoints object
#'
#' @param object Object created by call to [obsve()] or [matching_ve()]
#' @param alpha Significance level for confidence interval
#' @param seed Seed for reproducibility
#' @return A dataframe with lower and upper simultaneous confidence intervals for
#' each timepoint
#' @export
simultaneous_ci <- function(object, alpha, seed = NULL){

  stopifnot("Object must have more than 1 timepoint" = length(object$times) > 1)

  boot_samples <- object$boot_samples

  transformed_boot_samples <- list(logit(boot_samples[[1]]),
                                   logit(boot_samples[[2]]),
                                   log_ve(boot_samples[[3]]))
  names(transformed_boot_samples) <- names(boot_samples)

  #compute z_star for each term
  #set seed because call to get_z_star depends on random number generation
  if(!is.null(seed)){
    set.seed(seed)
  }

  z_star_results <- lapply(transformed_boot_samples, \(mat) get_z_star(mat, alpha = alpha))
  z_star <- sapply(z_star_results, \(x) x$z_star)
  excluded_timepoints <- lapply(z_star_results, \(x) as.numeric(x$excluded_timepoint))
  #print(z_star)

   #compute ci for each timepoint
   estimates <- object$estimates
   pt_est <- lapply(estimates, \(x) x[,"estimate"])


  risk_0_ci <- compute_wald_ci(x = pt_est[[1]],
                                     boot_x = boot_samples[[1]],
                                     transform = "logit",
                                     z_star = z_star[1])

  risk_1_ci <- compute_wald_ci(x = pt_est[[2]],
                                     boot_x = boot_samples[[2]],
                                     transform = "logit",
                                     z_star = z_star[2])

  ve_ci <- compute_wald_ci(x = pt_est[[3]],
                                 boot_x = boot_samples[[3]],
                                 transform = "log_ve",
                                 z_star = z_star[3])

  ci_estimates <- list(risk_0 = cbind(estimate = pt_est[[1]],  risk_0_ci),
                    risk_1 = cbind(estimate = pt_est[[2]],  risk_1_ci),
                        ve = cbind(estimate = pt_est[[3]],  ve_ci))

  #rename wald to simultaneous
  result <- lapply(seq_along(ci_estimates), \(i){
    x <- ci_estimates[[i]]
    # remove simultaneous CI estimates for excluded timepoints
    unused_times <- excluded_timepoints[[i]]
    x[rownames(x) %in% unused_times, ] <- NA
    colnames(x) <- gsub("wald", "simul", colnames(x))
    x[, c("estimate", "simul_lower", "simul_upper", "simul_n")]
  })
  names(result) <- names(ci_estimates)

  list(estimates = result,
       boot_samples = transformed_boot_samples,
       z_star = z_star,
       excluded_timepoints = excluded_timepoints)
}



#' Calculate critical value to use for simultaneous confidence bands
#'
#' If bootstrap matrix contains NA values at certain timepoints,
#' simultaneous confidence bands are only computed for timepoint regions where
#' at least 95% of the expected bootstrap samples are available
#'
#' @param mat A matrix of bootstrap estimates for a single term where rows
#' are bootstrap iterations and columns are time-points
#' @param alpha Significance level for confidence interval
#'
#' @return A list containing the critical value (scalar) and any timepoints
#' excluded from simultaneous confidence band region
#' @noRd
get_z_star <- function(mat, alpha ){
  #set infinite values to missing
  mat_clean <- mat
  mat_clean[is.infinite(mat)] <- NA

  #remove columns with a large proportion of  missing estimates
  bad_inds <- which(apply(mat_clean, 2, \(x) sum(is.na(x))/length(x)) > .05)


  if(length(bad_inds) > 0){
    times_removed <- colnames(mat_clean)[bad_inds]
    warning("Timepoints excluded: ", paste0(times_removed, collapse = ", "))
  }else{
    times_removed <- NULL
  }


  clean_inds <- setdiff(seq_len(ncol(mat_clean)), bad_inds)
  mat_sub <- mat_clean[,  clean_inds]
  covariance_sub <- stats::cov(mat_sub, use = "complete.obs")



  sds <- sqrt(diag(covariance_sub))
  estimate_draws <- MASS::mvrnorm(n = 10000,
                                  mu = rep(0, ncol(covariance_sub)),
                                  Sigma = (covariance_sub))
  max_z <- apply(estimate_draws, 1, \(x) max(abs(x/sds), na.rm = TRUE))
  z_star <- stats::quantile(max_z, 1-alpha)
  list(z_star = unname(z_star),
       excluded_timepoints = times_removed)
}
