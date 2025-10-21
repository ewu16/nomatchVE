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
                                keep_models = TRUE){


    matched_adata <- clean_matched_data(matched_data = matched_data,
                                        outcome_time = outcome_time,
                                        outcome_status = outcome_status,
                                        exposure = exposure,
                                        exposure_time = exposure_time,
                                        tau = tau)

   out <- compute_km_ve(adata = matched_adata,
                        adata_outcome_name = "match_T",
                        adata_event_name =  paste0("match_", outcome_status),
                        adata_trt_name = paste0("match_", exposure),
                        eval_times = eval_times,
                        keep_models = keep_models)
   return(out)


}
#' Get the marginal cumulative incidence in the treated and untreated
#' groups based on Kaplan Meier estimation.
#'
#' @inheritParams compute_marginal_ve
#' @inherit compute_marginal_ve return
#'
#' @keywords internal

compute_km_ve <- function(adata,
                          adata_outcome_name,
                          adata_event_name,
                          adata_trt_name,
                          eval_times,
                          keep_models = TRUE){

    outcome <- paste0("survival::Surv(", adata_outcome_name, ",", adata_event_name, ")")
    km_fit <- survival::survfit(stats::reformulate(adata_trt_name, response = outcome) , data = adata)
    surv_probs <- summary(km_fit, times = eval_times)
    cuminc_0 <- 1 - surv_probs$surv[surv_probs$strata == levels(surv_probs$strata)[1]]
    cuminc_1 <- 1 - surv_probs$surv[surv_probs$strata == levels(surv_probs$strata)[2]]
    rr <- cuminc_1/cuminc_0
    ve <- 1 - rr

    est  <- cbind(cuminc_0, cuminc_1, "risk_ratio" = rr, "vaccine_effectiveness" = ve)

    rownames(est) <- eval_times
    out <- list(pt_estimates = est)

    if(keep_models){
        out$models <- km_fit
    }
    return(out)
}





