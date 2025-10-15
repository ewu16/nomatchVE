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
#' @param gp_list List of two dataframes named `g_dist` and `p_dist` representing
#' the marginalizing distributions to use.
#'
#' @return A list containing the following:
#' \describe{
#'  \item{ci_estimates}{A list of matrices containing the lower and upper confidence interval bounds and
#'  bootstrap standard error for each term. When `confint_type = "wald"`, the bootstrapped standard errors
#'  on the transformed scale are also included.}
#'  \item{n_success_boot}{The number of bootstrap samples used to compute confidence interval}
#'  \item{boot_samples}{A matrix containing estimates from all bootstrap replications. Rows
#'  represent bootstrap iterations, columns the term estimated.}
#' }
#' @keywords internal
#'
estimate_nomatch_ci <- function(data,
                        outcome_name, event_name, trt_name, time_name,
                        adjust_vars, weighting,
                        times, censor_time, tau,
                        boot_formula_0 = NULL, boot_formula_1 = NULL,
                        matched_data =  NULL,
                        gp_list,
                        confint_type,
                        limit_type,
                        n_boot,
                        pt_est = NULL,
                        alpha = 0.05,
                        return_boot = TRUE,
                        n_cores = 1){

    one_boot_args <- list(data = data,
                          outcome_name = outcome_name,
                          event_name = event_name,
                          trt_name = trt_name,
                          time_name = time_name,
                          adjust_vars = adjust_vars ,
                          weighting = weighting,
                          times = times,
                          censor_time = censor_time,
                          tau = tau,
                          boot_formula_0 = boot_formula_0,
                          boot_formula_1 = boot_formula_1,
                          matched_data = matched_data,
                          gp_list = gp_list,
                          limit_type = limit_type)

   estimate_bootstrap_ci(one_boot_function_name = "one_boot_nomatch",
                       one_boot_args = one_boot_args,
                       confint_type = confint_type,
                       n_boot = n_boot,
                       pt_est = pt_est,
                       alpha = alpha,
                       return_boot = return_boot,
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
                        outcome_name, event_name, trt_name, time_name,
                        adjust_vars, weighting,
                        times, censor_time, tau,
                        boot_formula_0,
                        boot_formula_1,
                        matched_data,
                        gp_list,
                        limit_type){

    # --------------------------------------------------------------------------
    # 1. Create bootstrapped sample(s)
    # --------------------------------------------------------------------------
    #bootstrap original data
    boot_inds <- sample(1:nrow(data), replace = TRUE)
    boot_data <-  data[boot_inds,]


    #bootstrap matched adata if needed
    if(weighting == "matched" && limit_type == "limit"){
        #print("Bootstrap matched")
        boot_match_ids <- sample(unique(matched_data$pair_id), replace = TRUE)
        boot_match_inds <- as.vector(sapply(boot_match_ids, \(pair_id) which(matched_data$pair_id == pair_id)))
        boot_match_data <- matched_data[boot_match_inds,]
        boot_match_data$pair_id <- rep(1:length(boot_match_ids), each = 2)
    }else{
        boot_match_data <- NULL
    }

    # --------------------------------------------------------------------------
    # 2. Set marginalizing distribution based on limit type
    # --------------------------------------------------------------------------
    if(weighting == "custom"){
        limit_type <- "fixed"
    }
    if(limit_type == "fixed"){
        boot_weighting <- "custom"
        boot_custom_weights  <- gp_list
    }else if(limit_type == "limit"){
        boot_weighting <- weighting
        boot_custom_weights  <- NULL
    }

    # --------------------------------------------------------------------------
    # 3. Compute VE for bootstrapped data
    # --------------------------------------------------------------------------

    boot_ve <- get_one_nomatch_ve(data = boot_data ,
                          outcome_name = outcome_name,
                          event_name = event_name,
                          trt_name = trt_name,
                          time_name = time_name,
                          adjust_vars = adjust_vars ,
                          weighting = boot_weighting,
                          custom_weights = boot_custom_weights,
                          matched_dist_options = matched_dist(matched_data = boot_match_data),
                          times = times,
                          censor_time = censor_time,
                          tau = tau,
                          formula_0 = boot_formula_0,
                          formula_1 = boot_formula_1,
                          return_models = FALSE,
                          return_gp_list = FALSE,
                          return_matching = FALSE)

    boot_ve$estimates

}



