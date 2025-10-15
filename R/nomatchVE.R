#' Main function to estimate vaccine effectiveness (VE) without matching
#'
#' @description
#' `nomatchVE()` estimates VE over time without matching by using a G-computation
#' approach. It fits two conditional hazard models, one for each exposure type (e.g. vaccine vs no vaccine).
#' The models are used to predict time- and covariate- specific
#' cumulative incidences. These cumulative incidences are then marginalized over
#' the observed distributions of exposure uptake and baseline covariates among the exposed (e.g. vaccinated) group.
#'
#'
#' @param data A data frame with one row per individual containing
#'   the columns named in `outcome_name`, `event_name`, `trt_name`, `time_name`,
#'   and any variables listed in `adjust_vars`.
#' @param outcome_name Name of the possibly-censored time-to-event variable, measured from the chosen time origin.
#' @param event_name Name of the event indicator. The underlying column should be numeric
#'   (`1` = event, `0` = censored).
#' @param trt_name Name of the exposure indicator. The underlying column should be numeric
#'   (`1` = exposed during follow-up, `0` = never exposed during follow-up). Exposure indicator
#'   should indicate whether an individual is exposed while uncensored and before experiencing the endpoint of interest.
#' @param time_name Name of the time to exposure, measured from the chosen time origin; use `NA` if not exposed.
#' @param adjust_vars Character vector of baseline covariates to adjust for when fitting the hazard models.
#' @param tau Time after exposure (vaccination)  to exclude as part of an "immune build-up" period. Must be non-negative.
#'   Common values may include 0 days (no delay), 7 days (one week), 14 days (two weeks) and should
#'   match the biological understanding of when vaccine-induced immunity develops.
#' @param times Numeric vector of time after exposure at which to evaluate VE. All values must be greater than tau and
#'   should reflect meaningful clinical timepoints (e.g., 30, 60, 90 days).
#' @param censor_time Time after exposure at which exposed
#'   individuals are administratively censored during model fitting. Default:
#'   `max(times)`. This is used to limit borrowing of information from beyond the time period of interest.
#' @param weighting Character string specifying the type of marginalizing weights to use.
#'   Either:
#'   * `"observed"` (default): use the empirical distribution of time to exposure and covariates among the exposed as the
#'   marginalizing weights. This provides close alignment with the weights implicitly used in matching.
#'   * `"custom"`: use user-specified weights provided in the `custom_weights` argument
#' @param custom_weights a `list(g_dist, p_dist)` providing custom weights for marginalization.
#' Must have the following format:
#'   * `g_dist`: data frame with columns
#'      * `time_name` (time of exposure),
#'      *  all variables in `adjust_vars`
#'      * `group_name` (unique identifier for each covariate-group),
#'      * `prob` (probability of exposure at the given time within the covariate-group. Should sum to 1 within each covariate-group)
#'   * `p_dist`: data frame with columns
#'      *  all variables in `adjust_vars`
#'      * `prob` (probability of covariate-group. Should sum to 1 over all covariate groups.)
#' @param confint_type Method for constructing confidence intervals.
#'   One of `"wald"`, `"percentile"`, or `"both"`.
#'   - `"wald"`: Computes Wald-style intervals using bootstrap standard errors.
#'   - `"percentile"`: Computes percentile bootstrap intervals.
#'   - `"both"`: Computes and returns both sets of intervals.
#'   Default: `"wald"`. See **Confidence intervals** section.
#' @param n_boot Number of bootstrap replicates for confidence intervals.
#'   Recommended to use at least 1000 for publication-quality results. Use smaller values (e.g., 10-100) for initial exploration.
#'   Default: `0` (no bootstrapping).
#' @param alpha Significance level for confidence intervals (Confidence level = 100*(1-`alpha`)%). Default: `0.05`.
#' @param return_models Logical; return fitted conditional hazard models?
#'   Default: `TRUE`.
#' @param return_boot Logical; return bootstrap draws in `boot_samples`?
#'   Default: `TRUE`. Must be set to `TRUE` if user plans to use [add_simultaneous_ci()]
#'   to obtain simultaneous confidence intervals.
#' @param n_cores Integer; parallel cores for bootstrapping. Passed to
#'   `parallel::mclapply` as `mc.cores`. On Unix-like OS only; not available on
#'   Windows. Default: `1`.
#'
#'
#' @return An object of class `vefit` containing:
#' \describe{
#'   \item{estimates}{List of matrices: `cuminc_0`, `cuminc_1`, `ve`.
#'     Each matrix has rows corresponding to `times` and columns containing point estimates and confidence intervals.}
#'   \item{gp_list}{List with dataframes `g_dist`, `p_dist` used for marginalization.}
#'   \item{model_0}{Fitted hazard model under no exposure.}
#'   \item{model_1}{Fitted hazard model under exposure .}
#'   \item{n_success_boot}{Numeric vector with length equal to `length(times)`:
#'   indicates number of successful bootstrap replications per timepoint.}
#'   \item{boot_samples}{(If `return_boot = TRUE`) list of matrices with bootstrap
#'     estimates; rows = bootstrap iterations, cols = `times`.}
#' }
#'
#' The `vefit` object has methods for [print()], [summary()], and [plot()].
#' Use [add_simultaneous_ci()] to add simultaneous confidence intervals.
#'

#' @details
#' **Modeling.** Two Cox models are fit: one over the chosen time scale among the
#' not-yet-exposed, and one over time-since-exposure among the exposed who remain
#' at-risk `tau` days after exposure. For the exposed hazard model, exposure time is included as a natural cubic
#' spline with 4 degrees of freedom. Predictions are combined via G-computation and
#' marginalized.
#'
#'**Marginalizing weights.** When `weighting = observed`, the marginalizing weights
#' are the empirical distributions of exposure times and covariates among the exposed who remain
#' at-risk `tau` days after exposure. These weights are returned in the `vefit` object under `gp_list`.
#' They can also be obtained prior to the call to `nomatchVE()` by calling `get_observed_weights()`.
#'
#'
#' **Confidence intervals.** Wald CIs are constructed on transformed scales:
#'   logit for cumulative incidence; `log(1 - VE)` for VE, using bootstrap SEs,
#'   then back-transformed.
#'
#' **Parallelization.** Bootstraps can be parallelized on Unix via [parallel::mclapply()]
#' by providing `n_cores` argument.
#

#' @export
#'
#' @examples
#' # Fit vaccine effectiveness model using simulated data
#' data(simdata)
#'
#' fit <- nomatchVE(
#'   data = simdata,
#'   outcome_name = "Y",
#'   event_name = "event",
#'   trt_name = "V",
#'   time_name = "D_obs",
#'   adjust_vars = c("x1", "x2"),
#'   times = seq(30, 180, by = 30),
#'   tau = 14,
#'   n_boot = 5,
#'   n_cores = 2
#' )
#'
#' # View basic results
#' fit$estimates

nomatchVE <- function(data,
                  outcome_name,
                  event_name,
                  trt_name,
                  time_name,
                  adjust_vars,
                  tau,
                  times,
                  censor_time = max(times),
                  weighting = c("observed", "custom"),
                  custom_weights = NULL,
                  confint_type = c("wald", "percentile", "both"),
                  n_boot = 0,
                  alpha = 0.05,
                  return_models = TRUE,
                  return_boot = TRUE,
                  n_cores = 1
                  ){

     call <- match.call()
     weighting <- match.arg(weighting)
     confint_type <- match.arg(confint_type)

    # --------------------------------------------------------------------------
    # 0 - Prep
    # --------------------------------------------------------------------------

     # Check data/inputs
     validate_ve_inputs(data, outcome_name, event_name, trt_name, time_name, adjust_vars, times, tau, censor_time)

     if(weighting == "custom"){
         validate_marginalizing_weights(custom_weights, time_name, adjust_vars)
     }


     if (n_boot == 0) {
         if (confint_type != "none") {
             message("n_boot = 0, setting confint_type = 'none'")
         }
         confint_type <- "none"
     }
     if (confint_type == "none") {
         n_boot <- 0
     }

     # --------------------------------------------------------------------------
     # 1 - Get original estimate
     # --------------------------------------------------------------------------

     estimation_args <- list(data = data,
                             outcome_name = outcome_name,
                             event_name = event_name,
                             trt_name = trt_name,
                             time_name = time_name,
                             adjust_vars = adjust_vars,
                             times = times,
                             tau = tau,
                             censor_time = censor_time,
                             weighting = weighting,
                             custom_weights = custom_weights
                            )

     original <- do.call("get_one_nomatch_ve", estimation_args)



    # --------------------------------------------------------------------------
    # 2 - Get bootstrap CI
    # --------------------------------------------------------------------------
     estimate_nomatch_ci_args <- c(estimation_args[!names(estimation_args) %in% c("custom_weights")] ,
                           list(pt_est = original$estimates,
                                gp_list = original$gp_list,
                                confint_type = confint_type,
                                limit_type = "limit",
                                n_boot = n_boot,
                                alpha = alpha,
                                return_boot = return_boot,
                                n_cores = n_cores))


     boot_inference <- do.call("estimate_nomatch_ci", estimate_nomatch_ci_args)


     # --------------------------------------------------------------------------
     # 3 - Final result
     # --------------------------------------------------------------------------
     pt_est <- original$estimates
     ci_est <-boot_inference$ci_estimates
     estimates <- list(cuminc_0 = cbind(estimate = pt_est[,1], ci_est[[1]]),
                       cuminc_1 = cbind(estimate = pt_est[,2], ci_est[[2]]),
                       ve = cbind(estimate = pt_est[,3], ci_est[[3]]))


     # --------------------------------------------------------------------------
     # 4 - Return
     # --------------------------------------------------------------------------

     out <- list(
         estimates = estimates,
         gp_list = original$gp_list,
         model_0 = original$model_0,
         model_1 = original$model_1,
         outcome_name = outcome_name,
         event_name = event_name,
         trt_name = trt_name,
         time_name = time_name,
         adjust_vars = adjust_vars,
         times = times,
         censor_time = censor_time,
         tau = tau,
         confint_type = confint_type,
         n_boot = n_boot,
         n_success_boot = boot_inference$n_success_boot,
         boot_error_inds = boot_inference$error_inds,
         boot_na_list =boot_inference$boot_na_list,
         alpha = alpha,
         call = call)

    if(return_boot){
        out$boot_samples <- boot_inference$boot_samples
    }

     out$method <- "nomatchVE (G-computation)"
     class(out) <- "vefit"


    return(out)
}


