#' Compute VE point estimate in matched data set
#'
#' @description
#' First, creates the analysis matched data set based on `matched_data`.
#' Then computes Kaplan Meier estimates of cumulative incidence  on the analysis data set.
#'
#' @inheritParams matching_ve
#'
#' @return The object returned by [compute_km_ve()]
#' @keywords internal
#'
get_one_matching_ve <- function(matched_data,
                                outcome_time,
                                outcome_status,
                                exposure,
                                exposure_time,
                                tau,
                                eval_times,
                                keep_models = TRUE,
                                keep_adata = TRUE){


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

   check_pt_estimates(out$pt_estimates)

   #add some additional info
   if(keep_adata){
       out$matched_adata <- matched_adata
   }

   return(out)

}
#' Get the marginal cumulative incidence in the treated and untreated
#' groups based on Kaplan Meier estimation.
#'
#' @inheritParams matching_ve
#' @param adata A data frame that represents the analysis data set of a clinical trial.
#' @param adata_outcome_name Character string specifying the time to event variable in `adata`. The
#' time should be the time to event from vaccination/matched index date
#' @param adata_event_name Character string specifying the event variable in `adata`
#' @param adata_trt_name Character string specifying the treatment variable in `adata`
#'
#' @return A list containing the following:
#' \describe{
#' \item{estimates}{A matrix of estimates. The columns of the matrix are the cumulative
#'  incidence/VE terms and the rows are the requested time points for evaluation.}
#' \item{models}{If `keep_models = TRUE`, the models used to compute risk and VE}
#' }
#'
#' @keywords internal

compute_km_ve <- function(adata,
                          adata_outcome_name,
                          adata_event_name,
                          adata_trt_name,
                          eval_times,
                          keep_models = TRUE){

    # Build formula and fit KM
    outcome <- paste0("survival::Surv(", adata_outcome_name, ",", adata_event_name, ")")
    km_fit <- survival::survfit(stats::reformulate(adata_trt_name, response = outcome) , data = adata)

    # Summarize at requested times
    surv_probs <- summary(km_fit, times = eval_times)

    ## survival probabilities fo each group (assumed to have 2 groups_
    surv_0 <- surv_probs$surv[surv_probs$strata == levels(surv_probs$strata)[1]]
    surv_1 <- surv_probs$surv[surv_probs$strata == levels(surv_probs$strata)[2]]

    #Compute cumulative incidences and derived effect measures
    cuminc_0 <- 1 - surv_0
    cuminc_1 <- 1 - surv_1
    rr <- cuminc_1/cuminc_0
    ve <- 1 - rr

    est  <- cbind(
        cuminc_0 = cuminc_0,
        cuminc_1 = cuminc_1,
        risk_ratio = rr,
        vaccine_effectiveness = ve
    )

    rownames(est) <- eval_times

    out <- list(
        pt_estimates = est,
        models = if(keep_models) km_fit else NULL
    )
}





