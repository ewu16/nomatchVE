#' Internal function to compute VE point estimate
#'
#' @description This is an internal function that performs the actual VE computation.
#' It is called by [nomatchVE()],[nomatchVE_advanced()], and [one_boot_nomatch()]. Handles
#' different weighting schemes in a unified way. Users should typically call
#' these functions rather than calling this function directly.
#'
#' @inheritParams nomatchVE_advanced

#' @param return_matching Logical; return matched datasets? Default is
#'   TRUE if `weighting = "matched"`. When `weighting != "matched"`, this is
#'   automatically set to `FALSE`.
#'
##' @return List with components:
#' \describe{
#'   \item{estimates}{Matrix: `cuminc_0`, `cuminc_1`, `ve` (rows = timepoints)}
#'   \item{model_0, model_1}{Fitted Cox models (if `return_models = TRUE`)}
#'   \item{gp_list}{Marginalizing distributions (if `return_gp_list = TRUE`)}
#'   \item{matched_data, matched_adata}{Matched datasets (if `weighting = "matched"`)}
#' }
#'
#' @keywords internal
#' @export

get_one_nomatch_ve <- function(data,
                       outcome_name,
                       event_name,
                       trt_name,
                       time_name,
                       adjust_vars,
                       times,
                       censor_time,
                       tau,
                       weighting = c("observed", "custom", "matched"),
                       custom_weights = NULL,
                       matched_dist_options = NULL,
                       formula_0 = NULL,
                       formula_1 = NULL,
                       return_models = TRUE,
                       return_gp_list = TRUE,
                       return_matching = TRUE){

    # --------------------------------------------------------------------------
    # 0 - Check inputs/set defaults
    # --------------------------------------------------------------------------

    weighting <- match.arg(weighting)
    is_matching_marg <- (weighting == "matched")
    is_observed_marg <- (weighting == "observed")
    is_custom_marg <- (weighting == "custom")
    # if matching not applicable, do not return matching info
    if(!is_matching_marg){
        return_matching <- FALSE
    }

    # --------------------------------------------------------------------------
    # 1 - Fit survival models
    # --------------------------------------------------------------------------
    if(methods::is(formula_0, "coxph")){
        fit_0 <-  formula_0
    }else{
        fit_0 <- fit_model_0(data, outcome_name, event_name, trt_name, time_name,
                             adjust_vars, formula_0)
    }

    if(methods::is(formula_1, "coxph")){
        fit_1 <- formula_1
    }else{
        # Note: model_1 depends on t0 (censoring) and tau (subset)
        fit_1 <- fit_model_1(data, outcome_name, event_name, trt_name, time_name,
                             adjust_vars, censor_time, tau, formula_1)
    }
    # --------------------------------------------------------------------------
    # 2 - Get/fit marginalizing distributions
    # --------------------------------------------------------------------------
    if(is_matching_marg){
        if(is.null(matched_dist_options$matched_data)){
            print("Creating matched cohort")
            matched_cohort <- match_rolling_cohort(data = data,
                                                   outcome_name = outcome_name,
                                                   trt_name = trt_name,
                                                   time_name = time_name,
                                                   matching_vars = adjust_vars,
                                                   id_name = matched_dist_options$id_name,
                                                   replace = matched_dist_options$replace,
                                                   seed = matched_dist_options$seed)
            matched_dist_options$matched_data <- matched_cohort[[1]]
        }
        matched_adata <- clean_matched_data(matched_dist_options$matched_data,
                                             outcome_name = outcome_name,
                                             event_name = event_name,
                                             trt_name = trt_name,
                                             time_name = time_name,
                                             tau = tau)

        gp_list <- get_matching_gp(matched_adata = matched_adata,
                                   outcome_name = outcome_name,
                                   trt_name = trt_name,
                                   time_name = time_name,
                                   adjust_vars = adjust_vars,
                                   tau = tau)
    }else if(is_observed_marg){
        gp_list <- get_observed_gp(data = data,
                                   outcome_name = outcome_name,
                                   trt_name = trt_name,
                                   time_name = time_name,
                                   adjust_vars = adjust_vars,
                                   tau = tau)
    }else if(is_custom_marg){
        gp_list <- custom_weights
    }


    # --------------------------------------------------------------------------
    # 3 - Compute VE
    # --------------------------------------------------------------------------

    estimates <-compute_ve(fit_0, fit_1, time_name, times, tau, gp_list)


    # --------------------------------------------------------------------------
    # Return items
    # --------------------------------------------------------------------------
    out <- list(estimates = estimates)

    if(return_models){
        out$model_0 <- fit_0
        out$model_1 <- fit_1
    }
    if(return_gp_list){
        out$gp_list <- gp_list
    }
    if(return_matching){
        out$matched_data <-  matched_dist_options$matched_data
        out$matched_adata <- matched_adata
    }

    return(out)

}


