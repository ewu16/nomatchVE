#' Estimate bootstrapped confidence intervals for matching-based VE estimator
#'
#' @inheritParams matching_ve
#' @param pt_est A matrix of estimates. The columns of the matrix are the cumulative
#'  incidence/VE terms and the rows are the requested timepoints for evaluation.

#' @return A list containing the following:
#' \describe{
#'  \item{ci_estimates}{A list of matrices containing the lower and upper confidence interval bounds and
#'  bootstrap standard error for each term. When `ci_type = "wald"`, the bootstrapped standard errors
#'  on the transformed scale are also included.}
#'  \item{n_success_boot}{The number of bootstrap samples used to compute confidence interval}
#'  \item{boot_samples}{A matrix containing estimates from all bootstrap replications. Rows
#'  represent bootstrap iterations, columns the term estimated.}
#' }
#'
#' @keywords internal
#'
estimate_matching_ci <- function(matched_data,
                                 outcome_time = outcome_time,
                                 outcome_status,
                                 exposure,
                                 exposure_time,
                                 method = "km",
                                 adjust = NULL,
                                 eval_times,
                                 censor_time = NULL,
                                 tau,
                                 separate = FALSE,
                                 ci_type = "wald",
                                 limit_type = "fixed",
                                 data = NULL,
                                 id_name = "ID",
                                 matching_vars = NULL,
                                 replace  = FALSE,
                                 boot_reps = 0,
                                 pt_est = NULL,
                                 alpha = 0.05,
                                 keep_boot_samples,
                                 n_cores = 1){

    one_boot_args <- list(matched_data = matched_data,
                          outcome_time = outcome_time,
                          outcome_status = outcome_status,
                          exposure = exposure,
                          exposure_time = exposure_time,
                          method = method,
                          adjust = adjust,
                          eval_times = eval_times,
                          censor_time = censor_time,
                          tau = tau,
                          separate = separate,
                          limit_type = limit_type,
                          data = data,
                          id_name = id_name,
                          matching_vars = matching_vars,
                          replace = replace,
                          keep_boot_samples = keep_boot_samples)

    estimate_bootstrap_ci(one_boot_function_name = "one_boot_matching",
                        one_boot_args = one_boot_args,
                        ci_type = ci_type,
                        boot_reps = boot_reps,
                        pt_est = pt_est,
                        alpha = alpha,
                        keep_boot_samples = keep_boot_samples,
                        n_cores = n_cores)
}


#'  Compute one bootstrap replicate of matching-based VE point estimate
#'
#' @inheritParams estimate_matching_ci
#'
#' @return  A matrix of bootstrapped estimates where the the columns of the matrix are the cumulative
#'  incidence/VE terms and the rows are the requested time points for evaluation.
#'
#' @keywords internal
#'
one_boot_matching <- function(matched_data,
                              outcome_time,
                              outcome_status,
                              exposure,
                              exposure_time,
                              method,
                              adjust = NULL ,
                              eval_times,
                              censor_time = NULL,
                              tau,
                              separate = FALSE,
                              limit_type = "fixed",
                              data = NULL,
                              id_name = "ID",
                              matching_vars = NULL,
                              replace = FALSE,
                              keep_boot_samples = TRUE){

    # --------------------------------------------------------------------------
    # 1. Create bootstrapped sample(s)
    # --------------------------------------------------------------------------
    if(limit_type == "fixed"){
        #bootstrap from fixed matched cohort
        boot_matched_ids <- sample(unique(matched_data$match_id), replace = TRUE)
        boot_matched_inds <- as.vector(sapply(boot_matched_ids, \(match_id) which(matched_data$match_id == match_id)))
        boot_matched_data <- matched_data[boot_matched_inds,]
        boot_matched_data$match_id <- rep(1:length(boot_matched_ids), each = 2)
    }else if(limit_type == "limit"){
        #cat("Limit: \n")
        stopifnot("Need to provide original data for limit matching confidence intervals" =
                      !is.null(data))

        boot_inds <- sample(1:nrow(data), replace = TRUE)
        boot_data <-  data[boot_inds,]
        boot_id_name <- paste0("boot_",id_name)
        boot_data[[boot_id_name]] <- 1:nrow(boot_data)
        boot_matched_cohort <- match_rolling_cohort(data = boot_data,
                                                  outcome_time = outcome_time,
                                                  exposure = exposure,
                                                  exposure_time = exposure_time,
                                                  id_name =  boot_id_name,
                                                  matching_vars = matching_vars,
                                                  replace = replace)
        boot_matched_data <- boot_matched_cohort[[1]]

    }

    # --------------------------------------------------------------------------
    # 2. Compute VE for bootstrapped data
    # --------------------------------------------------------------------------
    boot_matching_ve <- get_one_matching_ve(matched_data = boot_matched_data,
                                            outcome_time = outcome_time,
                                            outcome_status = outcome_status,
                                            exposure = exposure,
                                            exposure_time = exposure_time,
                                            method = method,
                                            adjust = adjust,
                                            eval_times = eval_times,
                                            censor_time = censor_time,
                                            tau = tau,
                                            separate = separate,
                                            keep_models = FALSE)

    boot_matching_ve$estimates
}




