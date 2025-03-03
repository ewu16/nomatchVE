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
#' @export
#'
get_one_matching_ve <- function(matched_data,
                                outcome_name,
                                event_name,
                                trt_name,
                                time_name,
                                method,
                                adjust,
                                times,
                                censor_time,
                                tau,
                                pair_censoring = TRUE,
                                separate = TRUE,
                                return_models = TRUE){
    matched_adata <- clean_matched_data(matched_data = matched_data,
                                        outcome_name = outcome_name,
                                        event_name = event_name,
                                        trt_name = trt_name,
                                        time_name = time_name,
                                        tau = tau,
                                        pair_censoring = pair_censoring)

    estimates <- compute_marginal_ve(adata = matched_adata,
                                     adata_outcome_name = "T_d",
                                     adata_event_name = paste0(event_name, "_d"),
                                     adata_trt_name = paste0(trt_name, "_d"),
                                     method = method,
                                     adjust = adjust,
                                     times = times,
                                     censor_time = censor_time,
                                     separate = separate,
                                     return_models = return_models)
    estimates
}





