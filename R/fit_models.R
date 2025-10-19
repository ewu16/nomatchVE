#' Fit survival models to estimate exposure-specific hazards for G-computation
#' approach
#'
#'
#' @description `fit_model_0()` fits a Cox model to estimate risk for unexposed
#' individuals on the original time scale. Includes all individuals, censoring
#' exposed individuals at their time of exposure. By default, model is adjusted for by
#' `<adjust_vars>`, included as simple linear terms.
#'
#' `fit_model_1()` fits a Cox model to estimate risk for exposed individuals on
#' the time scale of time since exposure. Includes exposed individuals who
#' remain at risk `tau` days after exposure. Individuals are additionally
#' censored at `censor_time` days after exposure to avoid extrapolation beyond
#' the time period of interest. By default, model is adjusted for
#' `<adjust_vars>`, included as simple linear terms, and exposure time is included as a
#' natural cubic spline with 4 degrees of freedom.
#'
#' @inheritParams nomatchVE
#'
#' @param formula_0 Optional right hand side of the formula for model 0. By default, uses `adjust_vars`.
#'
#' @param formula_1 Optional right hand side of the formula for model 1. By default, uses `adjust_vars`
#'   plus natural spline of vaccination day (4 df). Default `NULL`
#'
#' @param censor_time Time after exposure at which exposed
#'   individuals are censored during model fitting to prevent extrapolation. By default,
#'   no censoring is applied.
#'
#' @return A fitted `coxph` object with additional component `$data` containing
#'   the analysis dataset used for fitting:
#' - For `fit_model_0()`: includes the survival tuple `(Y`, `event`)  and
#' covariates adjusted for in model, where `Y` is the time from time origin
#' to first of endpoint, censoring or exposure time (for exposed individuals).
#' - For `fit_model_1()`: includes the survival tuple `(T1`, `event`), `<time_name>`,
#' and covariates adjusted for in model, where `T1` is the time from exposure to
#' endpoint or censoring, with additional censoring by `censor_time`. Only includes
#' exposed individuals at risk `tau` days after exposure.
#' @export
#'
fit_model_0 <- function(data, outcome_name, event_name, trt_name, time_name,
                        adjust_vars, formula_0 = NULL){

    #Define formula
    if(is.null(formula_0)){
        formula_0 <- stats::reformulate(adjust_vars)
    }

    #Extract covariates
    covars <- stats::get_all_vars(formula_0, data)

    check_reserved_vars(covars, c("Y", "event"), "Model covariates")

    #Extract column values
    outcome <- data[[outcome_name]]
    event <- data[[event_name]]
    D_obs <- data[[time_name]]
    V <- data[[trt_name]]

    #Define survival outcomes
    # for model_0, censor exposed individuals at time of exposure
    # (provided they were at risk for the endpoint at time of exposure)
    exposed_at_risk <- V == 1 & D_obs < outcome
    Y <- ifelse(exposed_at_risk, D_obs, outcome)
    event <- ifelse(exposed_at_risk, 0, event)

    #Make full survival dataset
    data_0 <- cbind(Y, event, covars)

    #Fit model
    formula_0_with_response <- stats::update(formula_0, survival::Surv(Y,event) ~ .)
    fit_0 <- survival::coxph(formula_0_with_response, data_0,
                             model = TRUE)
    fit_0$data <- data_0
    fit_0
}

#' @rdname fit_model_0
#' @export
fit_model_1 <- function(data, outcome_name, event_name, trt_name, time_name,
                        adjust_vars, tau, censor_time = NULL, formula_1 = NULL){

    if(is.null(censor_time)){
        censor_time <- Inf #effectively prevents additional censoring
    }

    #Define formula
    if(is.null(formula_1)){
        d_term <- paste0("splines::ns(", time_name, ", df = 4)")
        formula_1 <-stats::reformulate(c(adjust_vars, d_term))
    }

    #Extract covariates
    covars <- stats::get_all_vars(formula_1, data)

    # Check for name conflicts for variables that will be added
    check_reserved_vars(covars, c("T1", "event"), "Model covariates")

    #Extract column values
    outcome <- data[[outcome_name]]
    event <- data[[event_name]]
    D_obs <- data[[time_name]]
    V <- data[[trt_name]]

    #Define survival outcomes
    #for model_1, time-to-event variable is time from vaccination
    T1 <- outcome - D_obs

    #for model_1, censor events after censor_time
    censor_flag <- T1 > censor_time
    T1 <- ifelse(censor_flag, censor_time, T1)
    event <- ifelse(censor_flag > censor_time, 0, event)

    #Create survival dataset
    data_full <- cbind(T1, event, covars)
    #for model_1, only include those vaccinated and at risk tau days after vaccination
    data_1 <- subset(data_full, V == 1 & T1 > tau)

    #Fit model
    formula_1_with_response <- stats::update(formula_1, survival::Surv(T1, event) ~ .)
    fit_1 <- survival::coxph(formula_1_with_response, data_1, model = TRUE)

    if(any(is.na(stats::coef(fit_1)))){
        warning("In model_1, at least one coefficient estimate is  NA",
                "Check model specification or sample size.")
    }
    fit_1$data <- data_1
    fit_1
}
