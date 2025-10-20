#' Compute a matching-based estimator of VE with confidence intervals
#'
#' @description This function is the main function for computing a matching-based estimator.
#'
#' @inheritParams nomatchVE
#' @inheritParams clean_matched_data
#' @param matched_data A data frame for the matched cohort

#' @return A list containing the following:
#' \describe{
#'  \item{estimates}{A list of matrices of the estimates at each timepoint. Rows of
#'  each matrix are the terms "cuminc_0", "cuminc_1", "vaccine_effectiveness". Columns of each matrix
#'  gives the point estimate and confidence intervals at the specified time point.}
#'  \item{eval_times}{The timepoints at which VE was evaluated}
#'  \item{n_success_boot}{A numeric vector of the number of successful bootstrap samples for each time point.(Success bootstrap samples are
#'  those that result in non-missing valid point estimates.)}
#'  \item{boot_samples}{If `keep_boot_samples = TRUE`, a list of matrices for each term that contain the bootstrap estimates where the rows are the bootstrap iterations and
#'  the columns are the time points.}
#' }
#' @export
#'
matching_ve <- function(matched_data,
                        outcome_time,
                        outcome_status,
                        exposure,
                        exposure_time,
                        tau,
                        eval_times,
                        effect = c("vaccine_effectiveness", "risk_ratio"),
                        ci_type = c("wald", "percentile", "both"),
                        boot_reps = 0,
                        alpha = 0.05,
                        keep_models = TRUE,
                        keep_boot_samples = TRUE,
                        n_cores = 1){

    call <- match.call()
    effect = match.arg(effect)
    ci_type <- match.arg(ci_type)


    # Check data/inputs
    stopifnot("<outcome_time> not in data" = outcome_time %in% names(matched_data))
    stopifnot("<outcome_status>  not in data" = outcome_status %in% names(matched_data))
    stopifnot("<exposure> not in data" = exposure %in% names(matched_data))
    stopifnot("<exposure_time> not in data" = exposure %in% names(matched_data))

    adjust <- NULL
    # --------------------------------------------------------------------------
    # 1 - Get original estimate
    # --------------------------------------------------------------------------
    estimation_args <- list(matched_data = matched_data,
                            outcome_time = outcome_time,
                            outcome_status = outcome_status,
                            exposure = exposure,
                            exposure_time = exposure_time,
                            tau = tau,
                            eval_times = eval_times)

    original <- do.call("get_one_matching_ve", estimation_args)


    # --------------------------------------------------------------------------
    # 2 - Get bootstrap CI
    # --------------------------------------------------------------------------
    estimate_matching_ci_args <- c(estimation_args,
                                   list(ci_type = ci_type,
                                        limit_type = "fixed",
                                        data = data,
                                        boot_reps = boot_reps,
                                        pt_est = original$estimates,
                                        alpha = alpha,
                                        keep_boot_samples = keep_boot_samples,
                                        n_cores = n_cores))

    boot_inference <- do.call("estimate_matching_ci", estimate_matching_ci_args)

    # --------------------------------------------------------------------------
    # 3 - Format and return results
    # --------------------------------------------------------------------------
    pt_est <- original$estimates
    ci_est <-boot_inference$ci_estimates

    # Only keep cumulative incidences and requested effect measure
    cuminc  <- list(cuminc_0 =  cbind(estimate = pt_est[, "cuminc_0"], ci_est[["cuminc_0"]]),
                    cuminc_1 = cbind(estimate = pt_est[, "cuminc_1"], ci_est[["cuminc_1"]]))

    effect_measure  <- setNames(list(cbind(estimate = pt_est[, effect], ci_est[[effect]])),
                        effect)

    estimates <- c(cuminc, effect_measure)
    boot_inference$boot_samples <-  boot_inference$boot_samples[names(estimates)]


    # --------------------------------------------------------------------------
    # 4 - Return
    # --------------------------------------------------------------------------

    out <- list(estimates = estimates,
                outcome_time = outcome_time,
                outcome_status = outcome_status,
                exposure = exposure,
                exposure_time = exposure_time,
                tau = tau,
                eval_times = eval_times,
                effect = effect,
                ci_type = ci_type,
                boot_reps = boot_reps,
                n_success_boot = boot_inference$n_success_boot,
                boot_error_list = boot_inference$error_list,
                boot_na_list =boot_inference$ boot_na_list,
                alpha = alpha,
                call = call)


    if(keep_models){
        out$models <- original$models
    }
    if(keep_boot_samples){
        out$boot_samples <- boot_inference$boot_samples
    }

    #for debugging
    #out$original <- original
    #out$one_boot_args <- boot_inference$one_boot_args

    out$method <- "matching"
    class(out) <- "vefit"


    return(out)
}




