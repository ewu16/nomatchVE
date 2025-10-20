#' Internal function to compute cumulative incidence and effect measures
#'
#' @description This is an internal function that performs the actual estimation
#' of cumulative incidences and derived effect measures
#' using the G-computation style approach. It is called by [nomatchVE()] and [one_boot_nomatch()].
#' Users should typically call these functions rather than calling this function directly.
#' For historical reasons, this function handles more complex inputs than
#' exposed in the `nomatchVE()` interface. In particular, it includes
#' an options to
#' - set censor time in hazard model for the exposed
#' - weight by the weights from a matched dataset and
#' - provide hazard model formulas or prefit objects
#'
#'
#' @inheritParams nomatchVE
#'
#' @param censor_time Time after exposure at which exposed
#'   individuals are administratively censored during model fitting. Default:
#'   `max(eval_times)`. This limits estimation to the observed  period of interest and
#'   prevents extrapolation beyond the largest evaluation time.
#'   Typically, users should leave this at the default.
#' @param weights_source Character string specifying how to construct marginalizing weights.
#'   One of:
#'   * `"observed"` (default): estimate weights from the observed data;
#'   * `"custom"`: use user-specified weights provided via `custom_weights` argument;
#' @param formula_0 Either a formula or a pre-fit `coxph` model object for the
#'   unvaccinated hazard (original time scale). If a fitted model is provided,
#'   it will be used directly instead of fitting a new model. Intended for advanced use.
#' @param formula_1 Either a formula or a pre-fit  `coxph` model object for the
#'   vaccinated hazard (time-since-exposure scale). If a fitted model is provided,
#'   it will be used directly instead of fitting a new model. Intended for advanced use.
#' @param return_gp_list Logical; return marginalizing weights? Default is
#'   TRUE.
#' @param return_matching Logical; return matched datasets? Default is
#'   TRUE if `weights_source = "matched"`. When `weights_source != "matched"`, this is
#'   automatically set to `FALSE`.
#'
##' @return List with components:
#' \describe{
#'   \item{estimates}{List of matrices: `cuminc_0`, `cuminc_1`, `risk_ratio`, `vaccine_effectiveness`
#'   where each matrix has one row per timepoint}
#'   \item{model_0, model_1}{Fitted Cox models (if `keep_models = TRUE`)}
#'   \item{gp_list}{Marginalizing distributions (if `return_gp_list = TRUE`)}
#'   \item{matched_data, matched_adata}{Matched datasets (if `weights_source = "matched"`)}
#' }
#'
#' @keywords internal
#' @export

get_one_nomatch_ve <- function(data,
                       outcome_time,
                       outcome_status,
                       exposure,
                       exposure_time,
                       covariates,
                       tau,
                       eval_times,
                       censor_time = max(eval_times),
                       effect = c("vaccine_effectiveness", "risk_ratio"),
                       weights_source = c("observed", "custom"),
                       custom_weights = NULL,
                       keep_models = TRUE,
                       return_gp_list = TRUE){

    # --------------------------------------------------------------------------
    # 0 - Check inputs/set defaults
    # --------------------------------------------------------------------------

    #Normalize user inputs
    effect <- match.arg(effect)
    weights_source <- match.arg(weights_source)


    # --------------------------------------------------------------------------
    # 1 - Fit survival models
    # --------------------------------------------------------------------------
    fit_0 <- fit_model_0(data, outcome_time, outcome_status, exposure, exposure_time,
                         covariates)

    # Note: model_1 depends on censor time and tau (subset)
    fit_1 <- fit_model_1(data, outcome_time, outcome_status, exposure, exposure_time,
                         covariates, tau, censor_time)

    # --------------------------------------------------------------------------
    # 2 - Get/fit marginalizing distributions depending on weights_source type
    # --------------------------------------------------------------------------
    if(identical(weights_source, "observed")){
        gp_list <- get_observed_gp(data = data,
                                   outcome_time = outcome_time,
                                   exposure = exposure,
                                   exposure_time = exposure_time,
                                   covariates = covariates,
                                   tau = tau)
    }else if(identical(weights_source, "custom")){
        gp_list <- canonicalize_weights(custom_weights)
    }

    # --------------------------------------------------------------------------
    # 3 - Compute VE
    # --------------------------------------------------------------------------
    cuminc <- compute_psi_bar_times(fit_0, fit_1, exposure_time, eval_times, tau, gp_list$g_weights, gp_list = gp_list)
    rr <-  cuminc[, "cuminc_1"]/cuminc[, "cuminc_0"]
    ve <- 1 - rr
    estimates <- cbind(cuminc, "risk_ratio" = rr, "vaccine_effectiveness" = ve)

    # --------------------------------------------------------------------------
    # Return items
    # --------------------------------------------------------------------------
    out <- list(estimates = estimates)

    if(keep_models){
        out$model_0 <- fit_0
        out$model_1 <- fit_1
    }
    if(return_gp_list){
        out$gp_list <- gp_list
    }


    return(out)

}


