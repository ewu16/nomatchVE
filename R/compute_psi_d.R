#' Compute day- and covariate- specific risks psi_v(d,x)
#'
#' @description Given the fitted survival models, compute estimates of psi_v(d,x)
#' for all d,x of interest.
#'
#' @inheritParams obsve
#' @param fit_0 A fitted survival model for the unvaccinated group
#' @param fit_1 A fitted survival model for the vaccinated group
#' @param gp_list List of marginalizing distributions. Used to determine which x, d
#' will be used for prediction
#'
#' @return A data frame containing predicted risk for each (d,x) of interest under
#' vaccine and no vaccine
#' @export
#'
compute_psi_d <- function(fit_0, fit_1, time_name, t0, tau, gp_list){

    pred_0 <- predict_from_model_0(fit_0, time_name, t0, tau, gp_list)
    pred_1 <- predict_from_model_1(fit_1, time_name, t0, tau, gp_list)
    psi_d <- merge(pred_0, pred_1)
    psi_d
}



#' Compute psi_0(d,x) and psi_1(d,x), respectively.
#'
#' @description These functions compute the conditional day- and covariate-specific cumulative
#' incidences.
#' `predict_from_model_0()` calculates `psi_0(d,x)` by using x as covariates in prediction and obtaining survival
#' at d + t0 and d + tau days
#' `predict_from_model_1()` calculates `psi_1(d,x)` by using x and d as covariates in prediction
#' and obtaining survival at t0 days
#'
#'
#' @inheritParams compute_psi_d
#'
#' @return A data frame of predictions
#' @export
#'
predict_from_model_0 <- function(fit_0, time_name, t0, tau, gp_list){

    # data_0_covars <- stats::get_all_vars(formula(delete.response(terms(fit_0))), fit_0$data)
    # unique_vals
    # new_data_0 <- expand.grid(data_0_covars)
    #pred_0 <- subset(gp_list$g_dist, select = -group_name)
    pred_0 <- gp_list$g_dist
    d_plus_tau <- d_plus_t0 <- pred_0
    d_plus_tau$Y <- pred_0[[time_name]] + tau
    d_plus_t0$Y <- pred_0[[time_name]] + t0

    d_plus_tau$event <- 0
    d_plus_t0$event <- 0

    pred_0$surv_d_plus_tau <-stats::predict(fit_0, d_plus_tau, type = "survival")
    pred_0$surv_d_plus_t0 <- stats::predict(fit_0, d_plus_t0,  type = "survival")
    pred_0$risk_0 <- 1 - pred_0$surv_d_plus_t0/pred_0$surv_d_plus_tau

    pred_0
}

#' @rdname predict_from_model_0
predict_from_model_1 <- function(fit_1, time_name, t0, tau, gp_list){
    #pred_1 <- subset(gp_list$g_dist, select = -group_name)
    pred_1 <- gp_list$g_dist
    time_tau <- time_t0 <- pred_1
    time_tau$T1_censored <-  tau
    time_t0$T1_censored <-  t0

    time_tau$event <- 0
    time_t0$event <- 0

    pred_1$surv_tau <- stats::predict(fit_1, time_tau, type = "survival")
    pred_1$surv_t0 <-  stats::predict(fit_1, time_t0,  type = "survival")
    pred_1$risk_1 <- 1 - pred_1$surv_t0/pred_1$surv_tau

    pred_1
}
