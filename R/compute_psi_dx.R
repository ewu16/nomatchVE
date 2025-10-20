#' Compute day- and covariate- specific cumulative incidences
#'
#' @description
#'  Wrapper that calls internal functions for predicting cumulative incidences from
#'  fitted hazard models. Returns the predicted exposure-specific cumulative
#'  incidences side by side for each time- and covariate- pair in `newdata`.
#'
#' @details
#' Definitions of the cumulative incidences returned:
#' \deqn{\psi_0(t_0; d,x) = 1 - S_0(d+t_0; x)\,/\,S_0(d+\tau; x)}
#' \deqn{\psi_1(t_0; d,x) = 1 - S_1(t_0; d,x)\,/\,S_1(\tau; d,x)}
#'
#' where \eqn{d} represents exposure time, \eqn{x} represents covariates, and
#' \eqn{S_v} represents the survival probability from hazard model for exposure \eqn{v}.
#'
#' @param fit_0 A fitted model returned from [fit_model_0()]
#' @param fit_1 A fitted model returned from [fit_model_1()]
#' @param exposure_time Name of the time-to-exposure variable  in `newdata`.
#'   Used to compute \eqn{\psi_0(t_0; d,x)} where \eqn{d + \tau} and \eqn{d + t_0} are needed.
#' @param tau Delay period
#' @param t0  Time since exposure at which to evaluate cumulative incidence.
#' @param newdata New data at which to do predictions.
#'
#' @return A data frame with one row per row of `newdata` with the predicted
#' cumulative incidences and component survival probabilities:
#' - `psi_0_dx`, `surv_0_d_plus_tau`, `surv_0_d_plus_t0`
#' - `psi_1_dx`, `surv_1_tau`, `surv_1_t0`

#' @seealso [predict_from_model_0()], [predict_from_model_1()]
#' @export
#'
compute_psi_dx_t0 <- function(fit_0, fit_1, exposure_time, t0, tau, newdata){
    #check new data argument
    pred_0 <- predict_from_model_0(fit_0, exposure_time, t0, tau, newdata)
    pred_1 <- predict_from_model_1(fit_1, exposure_time, t0, tau, newdata)
    psi_dx <- merge(pred_0, pred_1)
    psi_dx
}

#' Compute conditional cumulative incidences from fitted Cox models
#'
#' @description These functions compute the conditional day- and covariate-specific cumulative
#' incidence for each row of `newdata`.
#'
#' - `predict_from_model_0()` computes \eqn{\psi_0(t_0; d, x) = 1 - S_0(d+t_0 \mid x)/S_0(d+\tau \mid x)}
#'   by calling `predict(fit_0, newdata, type = "survival")` at eval_times
#' `d + t0` and `d + tau`.
#'
#' `predict_from_model_1()` computes \eqn{\psi_1(t_0; d, x) = 1 - S_1(t_0 \mid d, x)/S_1(\tau \mid d, x)}
#'  `by calling `predict(fit_1, newdata, type = "survival")` at eval_times
#'  `t0` and `tau`.
#'
#'
#' @return The `newdata` data frame with three additional columns:
#'  -  `predict_from_model_0()`: `surv_0_d_plus_tau`, `surv_0_d_plus_t0`, and `psi_0_dx`
#'  -  `predict_from_model_1()`: `surv_1_tau`, `surv_1_t0`, and `psi_1_dx`
#'

#' @keywords internal
#' @export
#'
#' @examples
#' # Fit hazard model under no vaccine
#' fit_0 <- fit_model_0(
#'   data = simdata,
#'   outcome_time = "Y",
#'   outcome_status = "event",
#'   exposure = "V",
#'   exposure_time = "D_obs",
#'   covariates = covariates
#'  )
#'
#' # Fit hazard model under vaccine
#' fit_1 <- fit_model_1(
#'   data = simdata,
#'   outcome_time = "Y",
#'   outcome_status = "event",
#'   exposure = "V",
#'   exposure_time = "D_obs",
#'   tau = tau
#'  )
#'
#'
#' # Define dataset for prediction,
#' # e.g. vaccinated indiviudals at risk tau days after vaccination
#' newdata <- simdata[simdata$V == 1 & (simdata$Y - simdata$D_obs) > 14,]
#'
#' # Predict from hazard model under no vaccine
#' predict_from_model_0(
#'     fit_0,
#'     exposure_time = "D_obs",
#'     t0 = 90,
#'     tau = 14,
#'     newdata = newdata
#' )
#'
#' # Predict from hazard model under vaccine
#' predict_from_model_1(
#'   fit_1,
#'   exposure_time = "D_obs",
#'   t0 = 90,
#'   tau = 14,
#'   newdata = newdata
#' )

predict_from_model_0 <- function(fit_0, exposure_time, t0, tau, newdata){
    surv_vars <- all.vars(formula(fit_0)[[2]])
    time_var <- surv_vars[1]
    event_var <- surv_vars[2]

    # Build newdata to predict survival at d+tau and d+t0
    newdata_d_plus_tau              <- newdata
    newdata_d_plus_tau[[time_var]]  <- newdata[[exposure_time]] + tau
    newdata_d_plus_tau[[event_var]] <- 0

    newdata_d_plus_t0               <- newdata
    newdata_d_plus_t0[[time_var]]   <- newdata[[exposure_time]] + t0
    newdata_d_plus_t0[[event_var]]  <- 0

    # Predict and store survival predictions
    pred_0 <- newdata
    pred_0$surv_0_d_plus_tau <-stats::predict(fit_0, newdata_d_plus_tau, type = "survival")
    pred_0$surv_0_d_plus_t0 <- stats::predict(fit_0, newdata_d_plus_t0,  type = "survival")
    pred_0$psi_0_dx <- 1 - pred_0$surv_0_d_plus_t0/pred_0$surv_0_d_plus_tau

    pred_0
}

#' @rdname predict_from_model_0
#'
predict_from_model_1 <- function(fit_1, exposure_time, t0, tau, newdata){
    surv_vars <- all.vars(formula(fit_1)[[2]])
    time_var <- surv_vars[1]
    event_var <- surv_vars[2]

    # Build newdata to predict survival at tau and t0
    newdata_tau              <- newdata
    newdata_tau[[time_var]]  <- tau
    newdata_tau[[event_var]] <- 0

    newdata_t0               <- newdata
    newdata_t0[[time_var]]   <- t0
    newdata_t0[[event_var]]  <- 0

    #This should always be 1 but sometimes returns NaN so just hard code
    surv_1_tau <- stats::predict(fit_1, newdata_tau, type = "survival")
    if(any(surv_1_tau[is.finite(surv_1_tau)] != 1)){
        stop("Survival at tau not 1")
    }

    # Predict survivals
    pred_1 <- newdata
    pred_1$surv_1_tau <- 1
    pred_1$surv_1_t0 <-  stats::predict(fit_1, newdata_t0,  type = "survival")
    pred_1$psi_1_dx <- 1 - pred_1$surv_1_t0/pred_1$surv_1_tau

    pred_1
}
