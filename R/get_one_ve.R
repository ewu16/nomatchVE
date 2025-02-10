#' Compute VE point estimate
#'
#' @description This function is called by [obsve] as well as by [one_boot_ve].
#'
#' @inheritParams obsve
#' @param return_models Logical: Should survival models be returned?
#' @param return_gp_list Logical: Should marginalizing distributions be returned?
#' @param return_matching Logical: Should matched datasets be returned? Default is
#' TRUE if marginalizing_dist = "matched" and FALSE otherwise.
#'
#' @return A list containing the following:
#' \describe{
#'  \item{estimates}{A named numeric vector containing point estimates for
#'  `psi_bar_0, psi_bar_1, ve`}
#'  \item{model_0, model_1}{If `return_models = TRUE`, the fitted survival models
#'   for unvaccinated and vaccinated.}
#'  \item{gp_list}{If `return_gp_list = TRUE`, the list of marginalizing distributions used.}
#'  \item{matched_data}{If `return_matching = TRUE`, the matched cohort}
#'  \item{matched_adata}{If `return_matching = TRUE`, the analysis-eligible
#'   matched cohort used to estimate `gp_list`}
#' }
#'
#' @export
#'
get_one_ve <- function( data,
                        outcome_name,
                        event_name,
                        trt_name,
                        time_name,
                        adjust_vars,
                        marginalizing_dist,
                        matched_dist_options,
                        t0,
                        censor_time,
                        tau,
                        formula_0,
                        formula_1,
                        return_models = TRUE,
                        return_gp_list = TRUE,
                        return_matching = NULL){

    # default return_matching is TRUE if marginalizing distribution is "matched data"
    # force to be false is marginalizing distribution is not "matched"
    if(!is.list(marginalizing_dist) && marginalizing_dist == "matched"){
        if(is.null(return_matching)){
            return_matching <- TRUE
        }
    }else{
        return_matching <-  FALSE
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
                             adjust_vars, t0, censor_time, tau, formula_1)
    }
    # --------------------------------------------------------------------------
    # 2 - Get/fit marginalizing distributions
    # --------------------------------------------------------------------------

    if(is.list(marginalizing_dist)){
        gp_list <- marginalizing_dist
    }else{
        #Prepare matched_adata if needed
        if(marginalizing_dist == "matched"){
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
                matched_dist_options$matched_data <-   matched_cohort[[1]]
            }
            matched_adata <- clean_matched_data(matched_dist_options$matched_data,
                                                outcome_name = outcome_name,
                                                event_name = event_name,
                                                trt_name = trt_name,
                                                time_name = time_name,
                                                tau = tau,
                                                pair_censoring = matched_dist_options$pair_censoring)
        }
        gp_list <- get_marginalizing_distributions(data, outcome_name, trt_name,
                                                   time_name, adjust_vars, tau,
                                                   marginalizing_dist_type = marginalizing_dist,
                                                   matched_adata = matched_adata)
    }

    # --------------------------------------------------------------------------
    # 3 - Compute psi_v(d,x)
    # --------------------------------------------------------------------------
    psi_d <- compute_psi_d(fit_0, fit_1, time_name, t0, tau, gp_list)

    # --------------------------------------------------------------------------
    # 4 - Marginalize
    # --------------------------------------------------------------------------
    psi_bar <- compute_psi_bar(psi_d, gp_list)
    psi_bar_ve <- add_ve(psi_bar)


    # --------------------------------------------------------------------------
    # Return items
    # --------------------------------------------------------------------------
    out <- list(estimates = psi_bar_ve)

    if(return_models){
        out$model_0 <- fit_0
        out$model_1 <- fit_1
    }
    if(return_gp_list){
        out$gp_list <- gp_list
    }
    if(return_matching){
        # Note: if objects are NULL then these items will not get added to list
        out$matched_data <-  matched_dist_options$matched_data
        out$matched_adata <- matched_adata
    }

    return(out)

}


