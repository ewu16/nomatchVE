#' Main function to estimate vaccine effectiveness using a G-computation
#' style estimator
#'
#' @description The primary function to estimate vaccine effectiveness.
#' @param data A data frame containing pertinent information for VE estimation.
#' @param outcome_name Character string specifying the name of the outcome
#'   variable in `data`. Outcome variable should be a numeric representing follow-up time relative to study
#'   start d0.
#' @param event_name Character string specifying the name of the event indicator
#'   variable in `data`. Event indicator should be numeric-valued: 1=event, 0=censored.
#' @param trt_name Character string specifying name of the vaccination indicator
#'   variable in `data`. Vaccination indicator is an indicator of ever receiving vaccine during the
#'   follow-up period and should be numeric-valued: 1=yes, 0 = no.
#' @param time_name Character string specifying name of the vaccination time
#'   variable in `data`. Vaccination time should be NA if never received vaccine
#'   during follow-up.
#' @param adjust_vars A character vector containing the names of variables in
#'   `data` to adjust for.
#' @param marginalizing_dist Character string describing the type of estimated
#'   day/covariate distributions to marginalize over. Values are "observed" or
#'   "matched". Alternatively, for pre-specified marginalizing distributions can
#'   be a list of 2 dataframes named `g_dist` and `p_dist`.
#' @param matched_dist_options If `marginalizing_dist == "matched"`, a list of
#'   parameters controlling process of creating matched cohort and matched analysis
#'   dataset from which marginalizing distributions will be estimated. List
#'   can be created by [matched_dist()].
#' @param times A numeric vector containing the times (relative to time of vaccination) at which to return
#'   the cumulative incidence estimates.
#' @param censor_time The time at which vaccinated individuals are censored during model fitting.
#'    By default, this is set to `max(times)`.
#' @param tau The time excluded after vaccination to allow building up of
#'   immunity
#' @param formula_0 A formula or pre-fit model for estimating hazards in
#'   unvaccinated group (recommended for use only by power users)
#' @param formula_1 A formula or pre-fit model for estimating hazards in
#'   vaccinated group. (recommended for use only by power users)
#' @param ci_type Character string indicating which type of confidence interval
#'   to return ("wald", "percentile", "both")
#' @param limit_type Character string indicating whether the marginalizing
#'   distributions of interest are the estimated distributions in the "limit"
#'   or "fixed" sense. For a prespecified marginalizing_dist, limit_type is
#'   always fixed.
#' @param n_boot Number of bootstrap replicated used to compute confidence
#'   intervals
#' @param alpha Significance level used to compute confidence intervals.
#'   Confidence intervals have nominal level `1 - alpha`.
#' @param return_models Logical value: Should fitted survival models be returned?
#' @param return_boot Logical value: Should bootstrap estimates be returned?
#' @param n_cores Number of cores to use to when running bootstrapping procedure. Passed as `mc.cores` argument to `parallel::mclapply`.
#' If `n_cores > 1`, this parallelizes the bootstrapping procedure on Unix-like systems (not available on Windows).
#'
#' @return A list. In addition to information related to the call/arguments to `obsve`, contains the following:
#' \describe{
#'  \item{estimates}{A list of matrices for cumulative incidences and VE.
#'  Each matrix contains the point estimate and confidence intervals for the specified term.}
#'  \item{gp_list}{A list containing the distributions g(d|x) and p(x) used for marginalization}
#'  \item{model_0}{The fit object used to predict risk for unvaccinated group}
#'  \item{model_1}{The fit object used to predict risk for vaccinated group}
 # \item{outcome_name}{Name of the outcome variable used}
 # \item{event_name}{Name of the event variable used}
 # \item{trt_name}{Name of the vaccine status variable used}
 # \item{time_name}{Name of the vaccination timing variable used}
 # \item{adjust_vars}{Name(s) of covariates adjusted for}
 # \item{t0}{The timepoint at which VE was computed and vaccinated individuals were censored
 # \item{tau}{The delay period for assessing vaccination}
 # \item{n_boot}{The number of bootstrap replicates requested}
#'  \item{n_success_boot}{A numeric vector of the number of successful bootstrap samples for each time point.(Success bootstrap samples are
#'  those that result in non-missing valid point estimates.)}
#'  \item{boot_samples}{If `return_boot = TRUE`, a list of matrices for each term that contain the bootstrap estimates where the rows are the bootstrap iterations and
#'  the columns are the time points.}
 # \item{alpha}{Significance level of confidence intervals returned}
 # \item{call}{The call to `obsve`}
#' }


#' @export
#'
# @examples




obsve <- function(data,
                  outcome_name,
                  event_name,
                  trt_name,
                  time_name,
                  adjust_vars,
                  marginalizing_dist,
                  matched_dist_options = NULL,
                  times,
                  censor_time = max(times),
                  tau = 14,
                  ci_type = "wald",
                  limit_type = "fixed",
                  n_boot = 0,
                  alpha = 0.05,
                  formula_0 = NULL,
                  formula_1 = NULL,
                  return_models = TRUE,
                  return_boot = TRUE,
                  n_cores = 1
                  ){

     call <- match.call()

    # --------------------------------------------------------------------------
    # 0 - Prep
    # --------------------------------------------------------------------------

     # Check data/inputs
     ## TODO

     if(!is.list(marginalizing_dist) && marginalizing_dist == "matched"){
        stopifnot("must provide matched_dist_options argument" = !is.null(matched_dist_options))
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
                             marginalizing_dist = marginalizing_dist,
                             matched_dist_options = matched_dist_options,
                             times = times,
                             censor_time = censor_time,
                             tau = tau,
                             formula_0 = formula_0,
                             formula_1 = formula_1)

     original <- do.call("get_one_ve", estimation_args)



    # --------------------------------------------------------------------------
    # 2 - Get bootstrap CI
    # --------------------------------------------------------------------------

     # for bootstrapping, cannot use prefit model, need to fit model in each sample
     if(methods::is(formula_0, "coxph")){
         boot_formula_0 <- stats::update(stats::formula(formula_0), NULL ~ .)
     }else{
         boot_formula_0 <- formula_0
     }

     if(methods::is(formula_1, "coxph")){
         boot_formula_1 <- stats::update(stats::formula(formula_1), NULL ~ .)
     }else{
         boot_formula_1 <- formula_1
     }


     unused_args <- names(estimation_args) %in% c("formula_0", "formula_1", "matched_dist_options")
     estimate_ci_args <- c(estimation_args[!unused_args],
                           list(boot_formula_0 = boot_formula_0,
                                boot_formula_1 = boot_formula_1,
                                matched_data = original$matched_data,
                                pt_est = original$estimates,
                                gp_list = original$gp_list,
                                ci_type = ci_type,
                                limit_type = limit_type,
                                n_boot = n_boot,
                                alpha = alpha,
                                return_boot = return_boot,
                                n_cores = n_cores))


     boot_inference <- do.call("estimate_ci", estimate_ci_args)


     # --------------------------------------------------------------------------
     # 3 - Final result
     # --------------------------------------------------------------------------
     pt_est <- original$estimates
     ci_est <-boot_inference$ci_estimates
     #ci_est <-boot_inference$ci_estimates
     estimates <- list(risk_0 = cbind(estimate = pt_est[,1], ci_est[[1]]),
                       risk_1 = cbind(estimate = pt_est[,2], ci_est[[2]]),
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
         ci_type = ci_type,
         limit_type = limit_type,
         n_boot = n_boot,
         n_success_boot = boot_inference$n_success_boot,
         boot_error_inds = boot_inference$error_inds,
         boot_na_list =boot_inference$ boot_na_list,
         alpha = alpha,
         call = call)

    if(!is.list(marginalizing_dist) && marginalizing_dist == "matched"){
        out$matched_data <- original$matched_data
        out$matched_adata <- original$matched_adata
        out$replace <- matched_dist_options$replace
        out$seed <- matched_dist_options$seed
        out$pair_censoring <- matched_dist_options$pair_censoring
    }
    if(return_boot){
        out$boot_samples <- boot_inference$boot_samples
    }

     #for debugging
     out$original <- original
     out$one_boot_args <- boot_inference$one_boot_args


    return(out)
}


