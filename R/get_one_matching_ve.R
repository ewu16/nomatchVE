#' Compute VE point estimate in matched data set
#'
#' @description
#' First, creates the analysis matched data set based on `matched_data`.
#' Then calls [compute_marginal_ve()] on the analysis data set, specifying
#' the appropriate outcome and treatment arguments.
#'
#' @inheritParams matching_ve
#'
#' @return The object returned by [compute_marginal_ve()]
#' @keywords internal
#'
get_one_matching_ve <- function(matched_data,
                                outcome_time,
                                outcome_status,
                                exposure,
                                exposure_time,
                                tau,
                                eval_times,
                                method = "km",
                                censor_time = NULL,
                                adjust = NULL,
                                separate = TRUE,
                                keep_models = TRUE){


    matched_adata <- clean_matched_data(matched_data = matched_data,
                                        outcome_time = outcome_time,
                                        outcome_status = outcome_status,
                                        exposure = exposure,
                                        exposure_time = exposure_time,
                                        tau = tau)

    estimates <- compute_marginal_ve(adata = matched_adata,
                                     adata_outcome_name = "match_T",
                                     adata_event_name =  paste0("match_", outcome_status),
                                     adata_trt_name = paste0("match_", exposure),
                                     method = method,
                                     adjust = adjust,
                                     eval_times = eval_times,
                                     censor_time = censor_time,
                                     separate = separate,
                                     keep_models = keep_models)
    estimates
}





