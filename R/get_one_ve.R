#' Compute VE point estimate
#'
#' @description This function computes point estimates of vaccine effectiveness
#' based on the proposed method. Internally called by [obsve] and [one_boot_ve].
#'
#' @inheritParams obsve
#' @param return_models Logical:  Should survival models be returned?
#' @param return_gp_list Logical: Should marginalizing distributions be returned?
#' @param return_matching Logical: Should matched datasets be returned? Default is
#' TRUE if marginalizing_dist = "matched" and FALSE otherwise.
#'
#' @return A list containing the following:
#' \describe{
#'  \item{estimates}{A matrix of estimates where the the columns of the matrix are the cumulative
#'  incidence/VE terms and the rows are the requested time points for evaluation.}
#'  \item{model_0, model_1}{If `return_models = TRUE`, the fitted survival models
#'   for unvaccinated and vaccinated.}
#'  \item{gp_list}{If `return_gp_list = TRUE`, the list of marginalizing distributions used.}
#'  \item{matched_data}{If `return_matching = TRUE`, the matched cohort}
#'  \item{matched_adata}{If `return_matching = TRUE`, the analysis-eligible
#'   matched cohort}
#' }
#'
#' @export
#'
get_one_ve <- function(data,
                       outcome_name,
                       event_name,
                       trt_name,
                       time_name,
                       adjust_vars,
                       marginalizing_dist,
                       matched_dist_options,
                       times,
                       censor_time,
                       tau,
                       formula_0,
                       formula_1,
                       return_models = TRUE,
                       return_gp_list = TRUE,
                       return_matching = TRUE){

    # --------------------------------------------------------------------------
    # 0 - Check inputs/set defaults
    # --------------------------------------------------------------------------

    # handle marginalizing_dist which can be list or string
    if(is.list(marginalizing_dist)){
        #TODO: check marginalizing_dist is proper gp_list?
        gp_list <- marginalizing_dist
        is_observed_marg <-is_matching_marg <- FALSE
    }else{
        stopifnot((marginalizing_dist %in% c("observed", "matched")))
        is_matching_marg <- (marginalizing_dist == "matched")
        is_observed_marg <- (marginalizing_dist == "observed")
    }

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
                                             tau = tau,
                                             pair_censoring = matched_dist_options$pair_censoring)

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


