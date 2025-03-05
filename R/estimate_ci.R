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
#'  bootstrap standard error for each term. When `ci_type = "wald"`, the bootstrapped standard errors
#'  on the transformed scale are also included.}
#'  \item{n_success_boot}{The number of bootstrap samples used to compute confidence interval}
#'  \item{boot_samples}{A matrix containing estimates from all bootstrap replications. Rows
#'  represent bootstrap iterations, columns the term estimated.}
#' }
#' @export
#'
estimate_ci <- function(data,
                        outcome_name, event_name, trt_name, time_name,
                        adjust_vars, marginalizing_dist,
                        times, censor_time, tau,
                        boot_formula_0, boot_formula_1,
                        matched_data,
                        gp_list,
                        ci_type,
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
                          marginalizing_dist = marginalizing_dist,
                          times = times,
                          censor_time = censor_time,
                          tau = tau,
                          boot_formula_0 = boot_formula_0,
                          boot_formula_1 = boot_formula_1,
                          matched_data = matched_data,
                          gp_list = gp_list,
                          limit_type = limit_type)

   estimate_general_ci(one_boot_function_name = "one_boot_ve",
                       one_boot_args = one_boot_args,
                       ci_type = ci_type,
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
#' @inheritParams estimate_ci
#'
#' @return  A matrix of bootstrapped estimates where the the columns of the matrix are the cumulative
#'  incidence/VE terms and the rows are the requested time points for evaluation.
#'
#' @export
#'
one_boot_ve <- function(data,
                        outcome_name, event_name, trt_name, time_name,
                        adjust_vars, marginalizing_dist,
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
    if(!is.list(marginalizing_dist) && marginalizing_dist == "matched" && limit_type == "limit"){
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
    if(is.list(marginalizing_dist)){
        limit_type <- "fixed"
        boot_marginalizing_dist <- marginalizing_dist
    }else if(limit_type == "fixed"){
        boot_marginalizing_dist <- gp_list
    }else if(limit_type == "limit"){
        boot_marginalizing_dist <- marginalizing_dist
    }else{
        stop("not a valid limit_type")
    }

    # --------------------------------------------------------------------------
    # 3. Compute VE for bootstrapped data
    # --------------------------------------------------------------------------

    boot_ve <- get_one_ve(data = boot_data ,
                          outcome_name = outcome_name,
                          event_name = event_name,
                          trt_name = trt_name,
                          time_name = time_name,
                          adjust_vars = adjust_vars ,
                          marginalizing_dist = boot_marginalizing_dist,
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



