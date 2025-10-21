#' Compute one bootstrap replicate for G-computation method
#'
#' @description This function computes the point estimate for
#' a bootstrap replicate.
#'
#' @inheritParams estimate_nomatch_ci
#'
#' @return  The `pt_estimates` component returned by [get_one_nomatch_ve()],
#'   for a bootstrap sample.
#'
#' @keywords internal
#'
one_boot_nomatch <- function(data,
                        outcome_time, outcome_status, exposure, exposure_time,
                        covariates,
                        eval_times, tau,
                        custom_gp_list){

    # --------------------------------------------------------------------------
    # 1. Create bootstrapped sample(s)
    # --------------------------------------------------------------------------
    #bootstrap original data
    boot_inds <- sample(1:nrow(data), replace = TRUE)
    boot_data <-  data[boot_inds,]

    # --------------------------------------------------------------------------
    # 2. Compute VE for bootstrapped data
    # --------------------------------------------------------------------------
    boot_ve <- get_one_nomatch_ve(data = boot_data,
                          outcome_time = outcome_time,
                          outcome_status = outcome_status,
                          exposure = exposure,
                          exposure_time = exposure_time,
                          covariates = covariates ,
                          custom_gp_list = custom_gp_list,
                          eval_times = eval_times,
                          tau = tau,
                          keep_models = FALSE,
                          return_gp_list = FALSE)

    boot_ve$pt_estimates

}



