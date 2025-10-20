#' Estimate bootstrapped confidence intervals for proposed VE estimator
#'
#' @inheritParams nomatchVE
#' @param boot_formula_0 A formula for estimating hazards in
#'   unvaccinated group for bootstrap samples (not a model)
#' @param boot_formula_1 A formula for estimating hazards in
#'   vaccinated group for bootstrap samples (not a model)
#' @param matched_data A data frame representing the matched cohort
#' @param pt_est A matrix of estimates. The columns of the matrix are the cumulative
#'  incidence/VE terms and the rows are the requested timepoints for evaluation.
#' @param gp_list List of two dataframes named `g_weights` and `p_weights` representing
#' the marginalizing distributions to use.
#'
#' @return A list containing the following:
#' \describe{
#'  \item{ci_estimates}{A list of matrices containing the lower and upper confidence interval bounds and
#'  bootstrap standard error for each term. When `ci_type = "wald"`, the bootstrapped standard errors
#'  on the transformed scale are also included.}
#'  \item{n_success_boot}{The number of bootstrap samples used to compute confidence interval}
#'  \item{boot_samples}{A matrix containing estimates from all bootstrap replications. Rows
#'  represent bootstrap iterations, columns the term estimated.}
#' }
#' @keywords internal
#'
estimate_nomatch_ci <- function(data,
                        outcome_time, outcome_status, exposure, exposure_time,
                        covariates, weights_source,
                        eval_times, censor_time, tau,
                        gp_list,
                        ci_type,
                        limit_type = "limit",
                        boot_reps,
                        pt_est = NULL,
                        alpha = 0.05,
                        keep_boot_samples = TRUE,
                        n_cores = 1){

    one_boot_args <- list(data = data,
                          outcome_time = outcome_time,
                          outcome_status = outcome_status,
                          exposure = exposure,
                          exposure_time = exposure_time,
                          covariates = covariates ,
                          weights_source = weights_source,
                          eval_times = eval_times,
                          censor_time = censor_time,
                          tau = tau,
                          gp_list = gp_list,
                          limit_type = limit_type)

   estimate_bootstrap_ci(one_boot_function_name = "one_boot_nomatch",
                       one_boot_args = one_boot_args,
                       ci_type = ci_type,
                       boot_reps = boot_reps,
                       pt_est = pt_est,
                       alpha = alpha,
                       keep_boot_samples = keep_boot_samples,
                       n_cores = n_cores)

}






#' Compute one bootstrap replicate of VE point estimate
#'
#' @description This function creates a bootstrapped sample and computes the
#' corresponding point estimate
#'
#' @inheritParams estimate_nomatch_ci
#'
#' @return  A matrix of bootstrapped estimates where the the columns of the matrix are the cumulative
#'  incidence/VE terms and the rows are the requested time points for evaluation.
#'
#' @keywords internal
#'
one_boot_nomatch <- function(data,
                        outcome_time, outcome_status, exposure, exposure_time,
                        covariates, weights_source,
                        eval_times, censor_time, tau,
                        gp_list,
                        limit_type){

    # --------------------------------------------------------------------------
    # 1. Create bootstrapped sample(s)
    # --------------------------------------------------------------------------
    #bootstrap original data
    boot_inds <- sample(1:nrow(data), replace = TRUE)
    boot_data <-  data[boot_inds,]


    # --------------------------------------------------------------------------
    # 2. Set marginalizing distribution based on limit type
    # --------------------------------------------------------------------------
    if(weights_source == "custom"){
        limit_type <- "fixed"
    }
    if(limit_type == "fixed"){
        boot_weighting <- "custom"
        boot_custom_weights  <- gp_list
    }else if(limit_type == "limit"){
        boot_weighting <- weights_source
        boot_custom_weights  <- NULL
    }

    # --------------------------------------------------------------------------
    # 3. Compute VE for bootstrapped data
    # --------------------------------------------------------------------------

    boot_ve <- get_one_nomatch_ve(data = boot_data ,
                          outcome_time = outcome_time,
                          outcome_status = outcome_status,
                          exposure = exposure,
                          exposure_time = exposure_time,
                          covariates = covariates ,
                          weights_source = boot_weighting,
                          custom_weights = boot_custom_weights,
                          eval_times = eval_times,
                          censor_time = censor_time,
                          tau = tau,
                          keep_models = FALSE,
                          return_gp_list = FALSE)

    boot_ve$estimates

}



