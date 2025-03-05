#' Compute a matching-based estimator of VE with confidence intervals
#'
#' @description This function is the main function for computing a matching-based estimator.
#'
#' @inheritParams obsve
#' @inheritParams clean_matched_data
#' @param matched_data A data frame for the matched cohort
#' @param method Character string specifying method for survival estimation ("cox" for
#' Cox proportional hazards regression model or "km" for Kaplan-Meier estimation)
#' @param adjust If `method = "cox"`, a formula or  character vector containing the names of covariates in
#' `matched_data` to adjust for.
#' @param separate  If `method = cox` and `separate = TRUE`, Cox regression models are fit in the treated and untreated groups separately.
#' @param limit_type The default value is `limit_type = "fixed"`, wherein bootstrap samples are formed from the
#' original matched data set. If `limit_type = "limit", the matching procedure is performed in each
#' bootstrap iteration.
#' @param data If `limit_type = "limit"`, the original data used to constructed the matched cohort
#' @param id_name  If `limit_type = "limit"`, character string representing the individual identifier variable  in `data`
#' @param matching_vars  If `limit_type = "limit"`, a character vector containing the names of variables in `data` to match on.
#' @param replace  If `limit_type = "limit"`, logical: Should matching be done with replacement?
#'
#' @return A list containing the following:
#' \describe{
#'  \item{estimates}{A list of matrices of the estimates at each timepoint. Rows of
#'  each matrix are the terms "risk_0", "risk_1", "ve". Columns of each matrix
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
                        method = "km",
                        adjust = NULL,
                        pair_censoring = TRUE,
                        separate = TRUE,
                        times,
                        censor_time,
                        tau,
                        ci_type = "wald",
                        n_boot = 0,
                        alpha = 0.05,
                        limit_type = "fixed",
                        data = NULL,
                        id_name = "ID",
                        matching_vars = NULL,
                        replace = FALSE,
                        return_models = TRUE,
                        return_boot = TRUE,
                        n_cores = 1){


    # Check data/inputs
    stopifnot("<outcome_name> not in data" = outcome_name %in% names(matched_data))
    stopifnot("<event_name>  not in data" = event_name %in% names(matched_data))
    stopifnot("<trt_name> not in data" = trt_name %in% names(matched_data))
    stopifnot("<time_name> not in data" = trt_name %in% names(matched_data))


    # --------------------------------------------------------------------------
    # 1 - Get original estimate
    # --------------------------------------------------------------------------
    estimation_args <- list(matched_data = matched_data,
                            outcome_name = outcome_name,
                            event_name = event_name,
                            trt_name = trt_name,
                            time_name = time_name,
                            method = method,
                            adjust = adjust,
                            times = times,
                            censor_time = censor_time,
                            tau = tau,
                            pair_censoring = pair_censoring,
                            separate = separate)

    original <- do.call("get_one_matching_ve", estimation_args)


    # --------------------------------------------------------------------------
    # 2 - Get bootstrap CI
    # --------------------------------------------------------------------------
    estimate_matching_ci_args <- c(estimation_args,
                                   list(ci_type = ci_type,
                                        limit_type = limit_type,
                                        data = data,
                                        id_name = id_name,
                                        matching_vars = matching_vars,
                                        replace = replace,
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
    estimates <- list(risk_0 = cbind(estimate = pt_est[,1], ci_est[[1]]),
                      risk_1 = cbind(estimate = pt_est[,2], ci_est[[2]]),
                      ve = cbind(estimate = pt_est[,3], ci_est[[3]]))


    out <- list(estimates = estimates,
                times = times,
                n_success_boot = boot_inference$n_success_boot,
                boot_error_list = boot_inference$error_list,
                boot_na_list =boot_inference$ boot_na_list)


    if(return_models){
        out$models <- original$models
    }
    if(return_boot){
        out$boot_samples <- boot_inference$boot_samples
    }

    #for debugging
    #out$original <- original
    #out$one_boot_args <- boot_inference$one_boot_args


    return(out)
}




