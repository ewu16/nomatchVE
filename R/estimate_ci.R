#' Estimate bootstrapped confidence intervals for VE
#'
#' @description This function computes Wald and percentile bootstrapped confidence
#' intervals for the proposed VE estimator.
#'
#' @inheritParams obsve
#' @param boot_formula_0 A formula for estimating hazards in
#'   unvaccinated group for bootstrap samples (not a model)
#' @param boot_formula_1 A formula for estimating hazards in
#'   vaccinated group for bootstrap samples (not a model)
#' @param matched_data A data frame representing the matched cohort
#' @param pt_est A named numeric vector containing point estimates for
#'  `psi_bar_0, psi_bar_1, ve`
#' @param gp_list List of two dataframes named `g_dist` and `p_dist` representing
#' the marginalizing distributions to use.

#' @return A list containing the following:
#' \describe{
#'  \item{ci_estimates}{A matrix containing the lower and upper confidence interval bounds
#'  and bootstrapped standard eror. When `ci_type = "wald"`, the bootstrapped standard errors
#'  on the transformed scale are also included.}
#'  \item{n_success_boot}{The number of bootstrap samples used to compute confidence interval}
#'  \item{boot_samples}{A matrix containing estimates from all bootstrap replications. Rows
#'  represent bootstrap iterations, columns the term estimated.}
#' }
#' @export
#'
#
estimate_ci <- function(data,
                        outcome_name, event_name, trt_name, time_name,
                        adjust_vars, marginalizing_dist,
                        t0, censor_time, tau,
                        boot_formula_0, boot_formula_1,
                        matched_data,
                        gp_list,
                        ci_type,
                        limit_type,
                        n_boot,
                        pt_est = NULL,
                        alpha = 0.05,
                        return_boot = TRUE){

    if(n_boot == 0){
        return(NULL)
    }

    # --------------------------------------------------------------------------
    # 1. Get bootstrap samples
    # --------------------------------------------------------------------------
    boot_samples <- run_bootstrap(data = data,
                                  outcome_name = outcome_name,
                                  event_name = event_name,
                                  trt_name = trt_name,
                                  time_name = time_name,
                                  adjust_vars = adjust_vars,
                                  marginalizing_dist = marginalizing_dist,
                                  t0 = t0,
                                  censor_time = censor_time,
                                  tau = tau,
                                  boot_formula_0 = boot_formula_0,
                                  boot_formula_1 = boot_formula_1,
                                  matched_data = matched_data,
                                  gp_list = gp_list,
                                  limit_type = limit_type,
                                  n_boot = n_boot)

    boot_samples_clean <- boot_samples[apply(boot_samples, 1, \(x) !any(is.na(x))),]

    # --------------------------------------------------------------------------
    # 2. Compute confidence intervals
    # --------------------------------------------------------------------------

    # Compute confidence interval
    if(ci_type == "wald"){
        ci <- compute_wald_ci(pt_est, boot_samples_clean, alpha)
    }else if(ci_type == "percentile"){
        ci <- compute_percentile_ci(boot_samples_clean, alpha)

    }else if(ci_type == "both"){
        ci <- cbind(compute_wald_ci(pt_est, boot_samples_clean, alpha),
                    compute_percentile_ci(boot_samples_clean, alpha))
    }

    # --------------------------------------------------------------------------
    # 3. Return confidence intervals and additional bootstrap information
    # --------------------------------------------------------------------------
    boot_sd <- apply(boot_samples_clean, 2, stats::sd)
    ci_estimates <- cbind(ci, boot_sd)
    n_success_boot <- nrow(boot_samples_clean)

    out <- list(ci_estimates = ci_estimates,
                n_success_boot = n_success_boot)

    if(return_boot){
        out$boot_samples <- boot_samples
    }

    return(out)

}


#' Compute bootstrapped Wald or percentile confidence intervals
#'
#' @description
#' `compute_wald_ci()`  computes Wald-style bootstrapped confidence intervals for
#' `psi_bar_0`, `psi_bar_1` and `ve`. Confidence intervals for `psi_bar_0, psi_bar_1`
#' are computed on the logodds scale (log(x/(1-x))) then transformed back to the
#' original scale. The confidence interval for `ve` is computed on the scale of `log(1-x)`
#' and then back-transformed.
#'
#' `compute_percentile_ci()` computes percentile bootstrapped confidence intervals.
#'
#' @inheritParams estimate_ci
#' @param boot_samples A matrix containing estimates from all bootstrap replications. Rows
#'  represent bootstrap iterations. Columns are named `psi_bar_0`, `psi_bar_1` and `ve`.
#' @param z_star Numeric vector with same length as `pt_est`. Specify the critical value to use when computing the confidence
#'  interval for each term in `pt_est`. Useful when computing simultaneous confidence intervals.
#'
#' @return Matrix containing the lower and upper bounds. For Wald confidence intervals,'
#' transformed standard errors are also included.
#' @export
#'
compute_wald_ci <- function(pt_est, boot_samples, alpha, z_star = NULL){

    if(!is.null(z_star)){
       z_crit <- z_star
    }else{
        z_crit <-  rep(stats::qnorm(alpha/2), length(pt_est))
    }

    psi_bar_0_ci <- compute_wald_logodds_ci(pt_est["psi_bar_0"], boot_samples[,"psi_bar_0"], z_crit[1])
    psi_bar_1_ci <- compute_wald_logodds_ci(pt_est["psi_bar_1"], boot_samples[,"psi_bar_1"], z_crit[2] )
    ve_ci <- compute_wald_log_ci(pt_est["ve"], boot_samples[,"ve"], z_crit[3])

    ci <- rbind(psi_bar_0_ci, psi_bar_1_ci, ve_ci)
    rownames(ci) <- names(pt_est)

    return(ci)

}


#' Compute transformed Wald intervals for cumulative incidence/VE
#'
#' @description
#' `compute_wald_logodds_ci`returns confidence interval based on transforming
#' Wald confidence intervals on the scale of `log(x/(1-x))`
#' `compute_wald_log_ci`returns confidence interval based on transforming
#' Wald confidence intervals on the scale of `log(1-x)`
#'
#'
#' @param x Point estimate
#' @param boot_x Numeric vector containing bootstrapped point estimates
#' @param z_crit_x Critical value to use to construct Wald-style confidence interval
#'
#' @return Numeric vector containing the lower and upper confidence interval, and
#' standard error used to confidence interval
#' @export
#'
compute_wald_logodds_ci <- function(x, boot_x, z_crit_x){
    x_logodds  <- log(x/(1-x))
    boot_logodds <- log(boot_x)/(1-boot_x)
    boot_sd <- stats::sd(boot_logodds)
    trans_ci <- stats::plogis(x_logodds + c(1, -1)*z_crit_x*boot_sd)
    wald <- c(trans_ci, boot_sd)
    names(wald) <- c("wald_lower", "wald_upper","wald_sd")
    wald
}

#' @rdname compute_wald_logodds_ci
compute_wald_log_ci <- function(x, boot_x, z_crit_x){
    x_log <- log(1 - x)
    boot_log <- log(1 - boot_x)
    boot_sd <- stats::sd(boot_log)
    trans_ci <- -1*(exp(x_log + c(-1, 1)*z_crit_x*boot_sd) - 1)
    wald <- c(trans_ci, boot_sd)
    names(wald) <- c("wald_lower", "wald_upper","wald_sd")
    wald
}






#' @rdname compute_wald_ci
compute_percentile_ci <- function(boot_samples, alpha){
    lower <- apply(boot_samples, 2, \(x) stats::quantile(x, alpha/2))
    upper <- apply(boot_samples, 2, \(x) stats::quantile(x, 1 - alpha/2))
    ci <- cbind(lower, upper)
    colnames(ci) <- c("percentile_lower", "percentile_upper")

    return(ci)
}
