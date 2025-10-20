#' Add simultaneous confidence intervals to vaccine effectiveness fit
#'
#' @description
#' Computes simultaneous confidence intervals, which maintain the specified coverage
#' level across all evaluation timepoints jointly. This is
#' useful for making inferences about the entire VE curve.
#'
#' @param object An object of class `vefit` created by [nomatchVE()] or [matching_ve()]. Must
#'   * contain evaluations at multiple timepoints (`length(object$eval_times)  > 0`),
#'   * contain bootstrap samples (`keep_boot_samples = TRUE` when fitting).
#' @param seed Integer seed for random number generation to ensure reproducible
#'   critical values for simultaneous confidence intervals. Default is `NULL` (no seed set).
#'
#' @return The input `vefit` object with the following modifications:
#'   \describe{
#'    \item{estimates}{Each matrix gets additional columns describing the simultaneous confidence interval bounds and construction:
#'    `simul_lower`, `simul_upper`, `simul_n`}
#'     \item{simul_z_star}{Critical values used for each term (cuminc_0, cuminc_1, ve)}
#'     \item{simul_excluded_timepoints}{Timepoints excluded from simultaneous bands
#'       due to insufficient bootstrap samples}
#'   }
#'
#' @details
#' Critical values are computed using the bootstrap covariance structure across
#' timepoints. If any timepoint has more than 5% missing bootstrap samples, it
#' is excluded from the simultaneous band and a warning is issued.
#'
#' @export
#'
#' @examples
#' # Fit model with bootstrap samples
#' fit <- nomatchVE(
#'   data = simdata,
#'   outcome_time = "Y",
#'   outcome_status = "event",
#'   exposure = "V",
#'   exposure_time = "D_obs",
#'   covariates = c("x1", "x2"),
#'   eval_times = seq(30, 180, by = 30),
#'   tau = 14,
#'   boot_reps = 100,
#'   keep_boot_samples = TRUE
#' )
#'
#' # Add simultaneous CIs
#' fit_simul <- add_simultaneous_ci(fit, seed = 123)
#'
#' # Look at results
#' fit_simul$estimates
#'
#' # Visualize
#' plot(fit_simul, ci_type = "simul")

add_simultaneous_ci <- function(object, seed = NULL){

  if (!is.vefit(object)) {
    stop("Object must be a vefit object", call. = FALSE)
  }

  if (length(object$eval_times) <= 1) {
    stop("Object must have more than 1 timepoint for simultaneous CIs",
         call. = FALSE)
  }

  if (is.null(object$boot_samples)) {
    stop("Object must contain bootstrap samples. ",
         "Rerun with keep_boot_samples = TRUE",
         call. = FALSE)
  }



  # Use alpha from the object
  alpha <- object$alpha

  boot_samples <- object$boot_samples

  transformed_boot_samples <- list(logit(boot_samples[[1]]),
                                   logit(boot_samples[[2]]),
                                   log_ve(boot_samples[[3]]))
  names(transformed_boot_samples) <- names(boot_samples)

  #compute z_star for each term
  z_star_results <- lapply(transformed_boot_samples, \(mat) get_z_star(mat, alpha = alpha, seed = seed))
  z_star <- sapply(z_star_results, \(x) x$z_star)
  excluded_timepoints <- lapply(z_star_results, \(x) as.numeric(x$excluded_timepoint))
  #print(z_star)

   #compute ci for each timepoint
   estimates <- object$estimates
   pt_est <- lapply(estimates, \(x) x[,"estimate"])


  cuminc_0_ci <- compute_wald_ci(x = pt_est[[1]],
                                     boot_x = boot_samples[[1]],
                                     transform = "logit",
                                     z_star = z_star[1])

  cuminc_1_ci <- compute_wald_ci(x = pt_est[[2]],
                                     boot_x = boot_samples[[2]],
                                     transform = "logit",
                                     z_star = z_star[2])

  ve_ci <- compute_wald_ci(x = pt_est[[3]],
                                 boot_x = boot_samples[[3]],
                                 transform = "log_ve",
                                 z_star = z_star[3])

  ci_estimates <- list(cuminc_0 = cbind(estimate = pt_est[[1]],  cuminc_0_ci),
                       cuminc_1 = cbind(estimate = pt_est[[2]],  cuminc_1_ci),
                       ve = cbind(estimate = pt_est[[3]],  ve_ci))

  # Format simultaneous CI results for each term
  simul_ci_list <- lapply(seq_along(ci_estimates), \(i){
    x <- ci_estimates[[i]]
    # remove simultaneous CI estimates for excluded timepoints
    unused_times <- excluded_timepoints[[i]]
    x[rownames(x) %in% unused_times, ] <- NA
    colnames(x) <- gsub("wald", "simul", colnames(x))
    x[, c("simul_lower", "simul_upper", "simul_n")]
  })
  names(simul_ci_list) <- names(ci_estimates)


  # Add simultaneous CIs to estimates
  object$estimates[[1]] <- cbind(object$estimates[[1]], simul_ci_list[[1]])
  object$estimates[[2]] <- cbind(object$estimates[[2]], simul_ci_list[[2]])
  object$estimates[[3]] <- cbind(object$estimates[[3]], simul_ci_list[[3]])


  object$simul_z_star <- z_star
  object$simul_excluded_timepoints <- excluded_timepoints

  return(object)
}


#' Calculate critical value for simultaneous confidence bands
#'
#' @description
#' Computes the critical value for simultaneous confidence intervals by simulating
#' from the bootstrap covariance structure (10,000 draws from multivariate normal).
#'
#' @param mat Matrix of bootstrap estimates (rows = iterations, columns = timepoints)
#' @param alpha Significance level (e.g., 0.05 for 95% confidence)
#' @param seed Integer seed for reproducibility. Default is `NULL`
#'
#' @return List containing:
#'   \describe{
#'     \item{z_star}{Critical value (scalar)}
#'     \item{excluded_timepoints}{Timepoints excluded due to >5% missing bootstrap estimates, or `NULL`}
#'   }
#'
#' @keywords internal
get_z_star <- function(mat, alpha, seed = NULL){

  #set seed for reproducibility as need to use multivariate normal random number generation
  if(!is.null(seed)){
    set.seed(seed)
  }

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
