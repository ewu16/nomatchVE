#' Fit survival models to estimate hazards and predict risk
#'
#' @inheritParams nomatchVE
#'
#' @description
#' `fit_model_0()` fits a Cox model to estimate risk for unvaccinated individuals
#'  * this model includes all individuals, censoring vaccinated individuals
#'  at their time of vaccination
#'  * TODO: need to censor individuals at d_max + t0? d_max not that well defined.
#'
#' `fit_model_1()` fits a Cox model to estimate risk for vaccinated individuals
#' * this model only includes vaccinated individuals at risk `tau` days after their vaccination
#' * furthermore, individuals are censored at `t0` days after vaccination
#'
#'
#' @return The fitted survival object
#' @export
#'
fit_model_0 <- function(data, outcome_name, event_name, trt_name, time_name,
                        adjust_vars, formula_0 = NULL){

    #Define formula
    if(is.null(formula_0)){
        formula_0 <- stats::reformulate(adjust_vars)
    }

    #Extract column values
    outcome <- data[[outcome_name]]
    event <- data[[event_name]]
    D_obs <- data[[time_name]]
    V <- data[[trt_name]]

    #Define survival outcomes
    #for model_0, censor vaccinated at time of vaccination
    Y <- ifelse(V == 1, D_obs, outcome)
    event <- ifelse(V == 1, 0, event)

    #Extract covariates
    covars <- stats::get_all_vars(formula_0, data)

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

fit_model_1 <- function(data, outcome_name, event_name, trt_name, time_name,
                        adjust_vars, censor_time, tau, formula_1 = NULL){

    #Define formula
    if(is.null(formula_1)){
        d_term <- paste0("splines::ns(", time_name, ", df = 4)")
        formula_1 <-stats::reformulate(c(adjust_vars, d_term))
    }

    #Extract column values
    outcome <- data[[outcome_name]]
    event <- data[[event_name]]
    D_obs <- data[[time_name]]
    V <- data[[trt_name]]

    #Define survival outcomes
    #for model_1, time-to-event variable is time from vaccination
    T1 <- outcome - D_obs
    #for model_1, censor events after censor_time
    T1_censored <- ifelse(T1 > censor_time, censor_time, T1)
    event <- ifelse(T1 > censor_time, 0, event)

    #Extract covariates
    covars <- stats::get_all_vars(formula_1, data)

    #Create survival dataset
    data_full <- cbind(T1_censored, event, covars)
    #for model_1, only include those vaccinated and at risk tau days after vaccination
    data_1 <- subset(data_full, V == 1 & T1_censored > tau)

    #Fit model
    formula_1_with_response <- stats::update(formula_1, survival::Surv(T1_censored,event) ~ .)
    fit_1 <- survival::coxph(formula_1_with_response, data_1, model = TRUE)

    if(any(is.na(stats::coef(fit_1)))){
        print("In model_1, at least one coefficient estimate is  NA")
        print(summary(fit_1))
    }
    fit_1$data <- data_1
    fit_1
}
