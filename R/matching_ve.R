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
#'  each matrix are the terms "cuminc_0", "cuminc_1", "ve". Columns of each matrix
#'  gives the point estimate and confidence intervals at the specified time point.}
#'  \item{times}{The timepoints at which VE was evaluated}
#'  \item{n_success_boot}{A numeric vector of the number of successful bootstrap samples for each time point.(Success bootstrap samples are
#'  those that result in non-missing valid point estimates.)}
#'  \item{boot_samples}{If `return_boot = TRUE`, a list of matrices for each term that contain the bootstrap estimates where the rows are the bootstrap iterations and
#'  the columns are the time points.}
#' }
#' @export
#'
matching_ve <- function(matched_data,
                        outcome_name,
                        event_name,
                        trt_name,
                        time_name,
                        tau,
                        times,
                        censor_time = max(times),
                        confint_type = c("wald", "percentile", "both"),
                        n_boot = 0,
                        alpha = 0.05,
                        return_models = TRUE,
                        return_boot = TRUE,
                        n_cores = 1){

    call <- match.call()

    # Check data/inputs
    stopifnot("<outcome_name> not in data" = outcome_name %in% names(matched_data))
    stopifnot("<event_name>  not in data" = event_name %in% names(matched_data))
    stopifnot("<trt_name> not in data" = trt_name %in% names(matched_data))
    stopifnot("<time_name> not in data" = trt_name %in% names(matched_data))

    adjust <- NULL
    # --------------------------------------------------------------------------
    # 1 - Get original estimate
    # --------------------------------------------------------------------------
    estimation_args <- list(matched_data = matched_data,
                            outcome_name = outcome_name,
                            event_name = event_name,
                            trt_name = trt_name,
                            time_name = time_name,
                            tau = tau,
                            times = times,
                            censor_time = censor_time)

    original <- do.call("get_one_matching_ve", estimation_args)


    # --------------------------------------------------------------------------
    # 2 - Get bootstrap CI
    # --------------------------------------------------------------------------
    estimate_matching_ci_args <- c(estimation_args,
                                   list(confint_type = confint_type,
                                        limit_type = "fixed",
                                        data = data,
                                        #id_name = id_name,
                                        #matching_vars = matching_vars,
                                        #replace = replace,
                                        n_boot = n_boot,
                                        pt_est = original$estimates,
                                        alpha = alpha,
                                        return_boot = return_boot,
                                        n_cores = n_cores))

    boot_inference <- do.call("estimate_matching_ci", estimate_matching_ci_args)

    # --------------------------------------------------------------------------
    # 3 - Format and return results
    # --------------------------------------------------------------------------
    pt_est <- original$estimates
    ci_est <-boot_inference$ci_estimates
    estimates <- list(cuminc_0 = cbind(estimate = pt_est[,1], ci_est[[1]]),
                      cuminc_1 = cbind(estimate = pt_est[,2], ci_est[[2]]),
                      ve = cbind(estimate = pt_est[,3], ci_est[[3]]))


    out <- list(estimates = estimates,
                outcome_name = outcome_name,
                event_name = event_name,
                trt_name = trt_name,
                time_name = time_name,
                adjust = adjust,
                times = times,
                censor_time = censor_time,
                tau = tau,
                confint_type = confint_type,
                n_boot = n_boot,
                n_success_boot = boot_inference$n_success_boot,
                boot_error_list = boot_inference$error_list,
                boot_na_list =boot_inference$ boot_na_list,
                alpha = alpha,
                call = call)


    if(return_models){
        out$models <- original$models
    }
    if(return_boot){
        out$boot_samples <- boot_inference$boot_samples
    }

    #for debugging
    #out$original <- original
    #out$one_boot_args <- boot_inference$one_boot_args

    out$method <- "matching"
    class(out) <- "vefit"


    return(out)
}




