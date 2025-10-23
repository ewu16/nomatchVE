#'Main function to estimate marginal cumulative incidences and derived effect
#'measures without matching
#'
#'@description `nomatchVE()` estimates marginal cumulative incidences under
#'  exposure and no exposure using a G-computation approach. The method fits two
#'  conditional hazard models- one for each exposure group- and uses
#'  these models to predict time- and covariate- specific cumulative incidences.
#'  The predictions are then marginalized to compute overall
#'  (marginal) cumulative incidences. The resulting cumulative incidences can be
#'   summarized as risk ratio (RR = 1 - risk_exposed/risk_unexposed) or
#'   vaccine effectiveness estimates (VE = 1 - RR).
#'
#'@param data A data frame with one row per individual containing the columns
#'  named in `outcome_time`, `outcome_status`, `exposure`, `exposure_time`, and any
#'  variables listed in `covariates`.
#'@param outcome_time Name of the time-to-event/censoring variable. Time should
#'  be measured from a given time origin (e.g. study start, enrollment, or age)
#'  for all individuals.
#'@param outcome_status Name of the event indicator. The underlying column should be
#'  numeric (`1` = event, `0` = censored).
#'@param exposure Name of the exposure indicator. The underlying column should
#'  be numeric (`1` = exposed during follow-up, `0` = never exposed during
#'  follow-up).
#'@param exposure_time Name of the time to exposure, measured from the chosen time
#'  origin; use `NA` if not exposed. Time must be measured in the same units
#'  (e.g. days) as that used for  `outcome_time`.
#'@param covariates Character vector of covariates to adjust for when fitting
#'  the hazard models. These covariates should include all known confounders of
#'  exposure and censoring measured at the chosen time origin.
#'@param tau Non-negative numeric value specifying the time after exposure that
#'  should be excluded from the risk evaluation period. This argument is
#'  primarily intended for vaccination exposures, where it is common to exclude
#'  the time after vaccination when immunity is still building. Time must be
#'  measured in the same units as that used for `outcome_time` and `exposure_time`
#'  and should reflect the biological understanding of when vaccine-induced
#'  immunity develops (usually 1-2 weeks). For non-vaccine exposures, `tau` can
#'  be set to 0 (no delay period).
#'@param eval_times Numeric vector specifying the timepoints at which to compute
#'  cumulative incidence and the derived effect measures. The timepoints should
#'  be expressed in terms of time since exposure. All values must be greater
#'  than `tau` and and should correspond to clinically meaningful follow-up
#'  durations, such as 30, 60, or 90 days after exposure. A fine grid of
#'  timepoints (e.g., `eval_times = (tau+1):100`) can be provided if cumulative
#'  incidence curves over time are desired.
#'@param effect Character. Type of effect measure to compute and return,
#'  based on the estimated cumulative incidences. Either
#'  `"vaccine_effectiveness"` (default) or `"risk_ratio"`.
#'@param weights_source Character string specifying the type of marginalizing weights
#'  to use. Either:
#'   - `"observed"` (default): set the marginalizing weights to the empirical
#'  distribution of exposure eval_times and covariates among the exposed. This
#'  provides close alignment with the weights implicitly used in matching.
#'   - `"custom"`: use the user-specified weights provided in the `custom_weights` argument.
#'@param custom_weights a `list(g_weights, p_weights)` providing weights for
#'  marginalizing the time- and covariate-specific cumulative incidences. Must
#'  have the following format:
#'   - `g_weights`: data frame with columns
#'      *  all variables in `covariates`
#'      * `exposure_time` (time of exposure),
#'      * `prob` (probability of exposure at the given time within the covariate-group;
#'      should sum to 1 within each covariate-group)
#'   - `p_weights`: data frame with columns
#'      *  all variables in `covariates`
#'      * `prob` (probability of covariate-group; should sum to 1 over all covariate groups.)
#'@param ci_type Method for constructing bootstrap confidence intervals. One of
#'  `"wald"`, `"percentile"`, or `"both"`.
#'   - `"wald"` (default): Computes Wald-style intervals using bootstrap standard errors.
#'   See **Confidence intervals** section for details.
#'   - `"percentile"`: Computes percentile bootstrap intervals.
#'   - `"both"`: Computes and returns both sets of intervals.
#'
#'@param boot_reps Number of bootstrap replicates for confidence intervals.
#'  Recommended to use at least 1000 for publication-quality results. Use
#'  smaller values (e.g., 10-100) for initial exploration. Default: `0` (no
#'  bootstrapping).
#'@param alpha Significance level for confidence intervals (Confidence level =
#'  100*(1-`alpha`)%). Default: `0.05`.
#'@param keep_models Logical; return the two fitted hazard models used to compute
#' cumulative incidences?
#'  Default: `TRUE`.
#'@param keep_boot_samples Logical; return bootstrap samples? Default:
#'  `TRUE`. Must be set to `TRUE` if user plans to use [add_simultaneous_ci()]
#'  to obtain simultaneous confidence intervals.
#'@param n_cores Integer; parallel cores for bootstrapping. Passed to
#'  `parallel::mclapply` as `mc.cores`. On Unix-like OS only; not available on
#'  Windows. Default: `1`.
#'
#'
#'@return An object of class `vefit` containing:
#' \describe{
#'   \item{estimates}{List of matrices:
#'   \describe{
#'      \item{`cuminc_0`}{ marginal cumulative incidence under no exposure}
#'      \item{`cuminc_1`}{ marginal cumulative incidence under exposure}
#'      \item{`<effect>`}{the selected effect measure}
#'   }
#'      Each matrix has one row per value in `eval_times` and columns including the
#'     point estimate (`estimate`) and, when requested, confidence limits of the form
#'     (`{wald/percentile}_lower`, `{wald/percentile}_upper`). }
#'   \item{weights}{List with dataframes `g_weights`, `p_weights` specifying
#'   the marginalizing weights used for averaging over exposure times and covariates.}
#'   \item{model_0}{Fitted hazard model for the unexposed group.
#'   See **Modeling** section for details.}
#'   \item{model_1}{Fitted hazard model for the exposed group.
#'   See **Modeling** section for details.}
#'   \item{n_success_boot}{Integer vector indicating the
#'   number of successful bootstrap replications per timepoint.}
#'   \item{boot_samples}{(If `keep_boot_samples = TRUE`) List of bootstrap draws
#'   (stored as matrices) for each
#'   returned quantity with names mirroring those in `estimates` (i.e. `cuminc_0`, `cuminc_1`, `<effect>`).
#'   Rows index bootstrap replicates and columns index `eval_times`.}
#' }
#'
#' The `vefit` object has methods for [print()], [summary()], and [plot()].
#' Use [add_simultaneous_ci()] to add simultaneous confidence intervals.
#'
#'
#'@details
#'
#' **Modeling.** Two Cox proportional hazards models are fit to estimate
#' exposure-specific cumulative incidences. The first models the hazard of the
#' outcome over the chosen time scale among individuals who have not yet been
#' exposed, with exposed individuals censored at their time of exposure. The
#' second models the hazard over time since exposure among individuals who were
#' exposed and remain at risk `tau` days after exposure. Both models adjust for
#' the specified covariates to help control for confounding. The second model
#' also flexibly adjusts for exposure time (as a natural cubic spline with 4
#' degrees of freedom) to capture time-varying background risk. Predicted risks
#' from both models are then marginalized over the specified covariate and
#' exposure-time distributions to obtain G-computation style cumulative
#' incidence estimates.
#'
#'
#'**Marginalizing weights.** When `weights_source = "observed"`, the marginalizing weights
#'are the empirical distributions of exposure eval_times and covariates among the
#'exposed who remain at-risk `tau` days after exposure. These weights are
#'returned in the `vefit` object under `gp_list`. They can also be obtained
#'prior to the call to `nomatchVE()` by calling `get_observed_weights()`.
#'
#'
#' **Confidence intervals.** Wald CIs are constructed on transformed scales:
#'\eqn{\text{logit}} for cumulative incidence; \eqn{\log{RR}} for risk ratios,
#'\eqn{\log{1 - VE}} for vaccine effectiveness, using bootstrap SEs. These are
#'then back-transformed to the original scale.
#'
#' **Parallelization.** Bootstraps can be parallelized on Unix via [parallel::mclapply()]
#'by providing `n_cores` argument.
#'
#'@export
#'
#' @examples
#' # Fit vaccine effectiveness model using simulated data
#' data(simdata)
#'
#' fit <- nomatchVE(
#'   data = simdata,
#'   outcome_time = "Y",
#'   outcome_status = "event",
#'   exposure = "V",
#'   exposure_time = "D_obs",
#'   covariates = c("x1", "x2"),
#'   eval_times = seq(30, 180, by = 30),
#'   tau = 14,
#'   boot_reps = 5,
#'   n_cores = 2
#' )
#'
#' # View basic results
#' fit$estimates

nomatchVE <- function(data,
                  outcome_time,
                  outcome_status,
                  exposure,
                  exposure_time,
                  covariates,
                  tau,
                  eval_times,
                  effect = c("vaccine_effectiveness", "risk_ratio"),
                  weights_source = c("observed", "custom"),
                  custom_weights = NULL,
                  ci_type = c("wald", "percentile", "both"),
                  boot_reps = 0,
                  alpha = 0.05,
                  keep_models = TRUE,
                  keep_boot_samples = TRUE,
                  n_cores = 1
                  ){

    # --------------------------------------------------------------------------
    # 0 - Prep
    # --------------------------------------------------------------------------

    # Normalize user choices
    call <- match.call()

    effect <- match.arg(effect)
    weights_source      <- match.arg(weights_source)
    ci_type   <- match.arg(ci_type)

    # Validate inputs
    validate_ve_inputs(
        data = data,
        outcome_time = outcome_time,
        outcome_status = outcome_status,
        exposure = exposure,
        exposure_time = exposure_time,
        covariates = covariates,
        tau = tau,
        eval_times = eval_times
    )

     if(identical(weights_source, "custom") ){
         validate_marginalizing_weights(
             custom_weights = custom_weights,
             exposure_time      = exposure_time,
             covariates    = covariates
         )

         # Format weights to improve efficiency of internal calls
         custom_gp_list <- canonicalize_weights(
             weights = weights,
             exposure_time  = exposure_time,
             covariates    = covariates
         )
     }else{
        custom_gp_list <- NULL
     }

     descrip <- get_basic_descriptives_nomatch(data,
                                       outcome_time = outcome_time,
                                       outcome_status = outcome_status,
                                       exposure = exposure,
                                       exposure_time = exposure_time,
                                       tau = tau)

     # --------------------------------------------------------------------------
     # 1 - Get original estimate
     # --------------------------------------------------------------------------
     estimation_args <- list(data = data,
                             outcome_time = outcome_time,
                             outcome_status = outcome_status,
                             exposure = exposure,
                             exposure_time = exposure_time,
                             covariates = covariates,
                             tau = tau,
                             eval_times = eval_times,
                             custom_gp_list = custom_gp_list
                            )

     original <- do.call(get_one_nomatch_ve, estimation_args)
     pt_est <- original$pt_estimates

    # --------------------------------------------------------------------------
    # 2 - Get bootstrap CI
    # --------------------------------------------------------------------------
    # Helper returns NULL if boot_reps = 0
     boot_inference <- estimate_bootstrap_ci(
         one_boot_function  = one_boot_nomatch,
         one_boot_args      = estimation_args,
         ci_type            = ci_type,
         boot_reps          = boot_reps,
         pt_est             = pt_est,
         alpha              = alpha,
         keep_boot_samples  = keep_boot_samples,
         n_cores            = n_cores
     )

     ci_est <-boot_inference$ci_estimates

     # --------------------------------------------------------------------------
     # 3 - Combine estimates with bootstrap CI
     # --------------------------------------------------------------------------
     terms_keep <- c("cuminc_0", "cuminc_1", effect)

     add_ci_columns <- function(term, pt_est, ci_est){
         x <- cbind(estimate = pt_est[, term], ci_est[[term]])
         rownames(x) <- rownames(pt_est)
         x
     }
     estimates <- stats::setNames(
         lapply(terms_keep, \(term) add_ci_columns(term, pt_est, ci_est)),
         terms_keep
         )

     # --------------------------------------------------------------------------
     # 4 - Return
     # --------------------------------------------------------------------------

     # Define the weights to return
     if(identical(weights_source, "custom")){
         weights <- custom_weights
     }else{
         weights <- gp_to_weights(original$gp_list)
     }

    # Build return object
     out <- list(
         # Core output
         estimates = estimates,
         model_0 = original$model_0,
         model_1 = original$model_1,
         weights = weights,

         # Bootstrap information if available
         n_success_boot   = boot_inference$n_success_boot,
         boot_errors  = boot_inference$boot_errors,
         boot_nas     = boot_inference$boot_nas,
         boot_samples     = if (keep_boot_samples) boot_inference$boot_samples[terms_keep] else NULL,

         # User-provided or default specifications
         outcome_time = outcome_time,
         outcome_status = outcome_status,
         exposure = exposure,
         exposure_time = exposure_time,
         covariates = covariates,
         tau = tau,
         eval_times = eval_times,
         effect = effect,
         ci_type = ci_type,
         boot_reps = boot_reps,
         alpha = alpha,

         # Meta information
         call = call,
         descrip = descrip,
         method =  "nomatchVE (G-computation)"
         )

     class(out) <- "vefit"

    return(out)
}


