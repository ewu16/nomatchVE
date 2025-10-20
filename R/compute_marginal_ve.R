#' Get the marginal cumulative incidence in the treated and untreated
#' groups based on a Cox model or Kaplan Meier curve.
#'
#'#' @description
#' Internally calls [compute_cox_ve()] or [compute_km_ve()] depending
#' on the value of `method`
#'
#'
#'
#' @inheritParams matching_ve
#' @param adata A data frame that represents the analysis data set of a clinical trial.
#' @param adata_outcome_name Character string specifying the time to event variable in `adata`. The
#' time should be the time to event from vaccination/matched index date
#' @param adata_event_name Character string specifying the event variable in `adata`
#' @param adata_trt_name Character string specifying the treatment variable in `adata`
#'
#' @return A list containing the following:
#' \describe{
#' \item{estimates}{A matrix of estimates. The columns of the matrix are the cumulative
#'  incidence/VE terms and the rows are the requested time points for evaluation.}
#' \item{models}{If `keep_models = TRUE`, the models used to compute risk and VE}
#' }
#' @keywords internal
#'
compute_marginal_ve <- function(adata,
                                adata_outcome_name,
                                adata_event_name,
                                adata_trt_name,
                                eval_times,
                                method,
                                adjust = NULL,
                                censor_time = max(eval_times),
                                separate = TRUE,
                                keep_models = TRUE){
    stopifnot(method %in% c("km", "cox"))
    if(method == "km"){
        marginal_risk <- compute_km_ve(adata,
                                         adata_outcome_name,
                                         adata_event_name,
                                         adata_trt_name,
                                         eval_times,
                                         keep_models = keep_models)
    }else if(method == "cox"){
        marginal_risk <- compute_cox_ve(adata,
                                          adata_outcome_name,
                                          adata_event_name,
                                          adata_trt_name,
                                          adjust,
                                          eval_times,
                                          censor_time,
                                          separate = separate,
                                          keep_models = keep_models)
    }
    marginal_risk
}

#' Get the marginal cumulative incidence in the treated and untreated
#' groups based on Kaplan Meier estimation.
#'
#' @inheritParams compute_marginal_ve
#'
#' @return A list containing the following:
#' \describe{
#' \item{estimates}{A matrix of estimates. The columns of the matrix are the cumulative
#'  incidence/VE terms and the rows are the requested time points for evaluation.}
#' \item{models}{If `keep_models = TRUE`, the models used to compute risk and VE}
#'}
#'
#' @keywords internal

compute_km_ve <- function(adata,
                          adata_outcome_name,
                          adata_event_name,
                          adata_trt_name,
                          eval_times,
                          keep_models = TRUE){
    outcome <- paste0("survival::Surv(", adata_outcome_name, ",", adata_event_name, ")")
    km_fit <- survival::survfit(stats::reformulate(adata_trt_name, response = outcome) , data = adata)
    surv_probs <- summary(km_fit, times = eval_times)
    cuminc_0 <- 1 - surv_probs$surv[surv_probs$strata == levels(surv_probs$strata)[1]]
    cuminc_1 <- 1 - surv_probs$surv[surv_probs$strata == levels(surv_probs$strata)[2]]
    rr <- cuminc_1/cuminc_0
    ve <- 1 - rr

    estimates <- cbind(cuminc_0, cuminc_1,
                       "risk_ratio" = rr,
                       "vaccine_effectiveness" = ve)

    rownames(estimates) <- eval_times


    out <- list(estimates = estimates)

    if(keep_models){
        out$models <- km_fit
    }
    return(out)
}


#' Get the marginal cumulative incidence in the treated and untreated
#' groups based on Cox model(s)
#'
#' @description This function uses a Cox model to predict the marginal cumulative
#' incidence of an endpoint `t0` days after vaccination. The marginal
#' cumulative incidence is computed by predicting individual survival probabilities
#' and averaging these predictions over everyone in the treatment group.
#'
#' @inheritParams compute_marginal_ve
#'
#' @return A list containing the following:
#' \describe{
#' \item{estimates}{A matrix of estimates. The columns of the matrix are the cumulative
#'  incidence/VE terms and the rows are the requested time points for evaluation.}
#' \item{models}{If `keep_models = TRUE`, the models used to compute risk and VE}
#' }
#'
#' @keywords internal
#'
compute_cox_ve <- function(adata,
                             adata_outcome_name,
                             adata_event_name,
                             adata_trt_name,
                             adjust,
                             eval_times,
                             censor_time,
                             separate = TRUE,
                             keep_models = TRUE){
    # --------------------------------------------------------------------------
    # 0. Prep data
    # --------------------------------------------------------------------------

    #Define survival outcomes
    tt <- adata[[adata_outcome_name]]
    event <- adata[[adata_event_name]]
    trt <- adata[[adata_trt_name]]


    #Censor at censor_time if given
    if(!is.null(censor_time)){
        #use flag due to redefining tt which change logic for second if else
        flag <- tt > censor_time
        tt <- ifelse(flag, censor_time, tt)
        event <- ifelse(flag , 0, event)
    }

    #Get covariates and formulas for models
    if(methods::is(adjust, "formula")){
        covars <- stats::get_all_vars(adjust, adata)
        formula_RHS <- stats::reformulate(c(attr(stats::terms(adjust), "term.labels"), "trt"))

    }else{
        covars <- adata[,adjust]
        formula_RHS <- stats::reformulate(c(adjust, "trt"))
    }


    tmp_data <- cbind(tt, event, trt, covars, dummy = 1)

    # --------------------------------------------------------------------------
    # 1. Fit models
    # --------------------------------------------------------------------------

    if(separate == TRUE){
        split_var <- "trt"
        formula_RHS <- stats::update(formula_RHS,  ~ . - trt)
    }else{
        #dummy splitting factor
        split_var <- "dummy"
    }

    data_list <- split(tmp_data, tmp_data[[split_var]])

    model_formula <- stats::update(formula_RHS, survival::Surv(tt,event) ~ .)
    models <- lapply(data_list, \(df){
        model <- survival::coxph(model_formula, data = df, model = TRUE)
        if(model$nevent == 0){
            #lack of formula in object returned can cause issues with prediction
            model$formula <- model_formula
            print("In Cox regression model, no events were observed")
        }else if(any(is.na(stats::coef(model)))){
            print("In Cox regression model, at least one coefficient estimate is  NA")
            print(summary(model))
        }
        model
    })


    # --------------------------------------------------------------------------
    # 2. Make predictions
    # --------------------------------------------------------------------------

    # Do prediction and compute marginal risk
    predict_at_t <- function(model_data, model, t){
        pred_data <- model_data
        pred_data$tt <- t
        if(model$nevent == 0){
            #hard code predicted survival probability when no events (predict returns NA)
            surv_prob <- rep(1,nrow(pred_data))
        }else{
            surv_prob <- stats::predict(model, newdata = pred_data, type = "survival")
        }
        data.frame(t0 = t, trt = model_data$trt, surv_prob)
    }

    pred_at_times <- function(model_data, model, eval_times){
        pred_list <- lapply(eval_times, \(t) predict_at_t(model_data, model, t))
        names(pred_list) <- eval_times
        do.call("rbind", pred_list)
        rownames(pred_list) <- NULL
        pred_list
    }


    combine_list <- lapply(seq_along(models), \(i){
        pred_at_times(model_data = data_list[[i]],
                      model = models[[i]],
                      eval_times)
    })

    # --------------------------------------------------------------------------
    # 3. Marginalize
    # --------------------------------------------------------------------------


    full_pred <- do.call("rbind", unlist(combine_list, recursive = FALSE))
    rownames(full_pred) <- NULL

    #by using aggregate can keep things simpler for separate and together case
    cuminc_df <- stats::aggregate(surv_prob ~trt + t0, data = full_pred,
                                FUN =  \(x){1 - mean(x)},
                                na.action = stats::na.pass)


    # --------------------------------------------------------------------------
    # 4. Format result
    # --------------------------------------------------------------------------
    cuminc_0 <- cuminc_df$surv_prob[cuminc_df$trt == 0]
    cuminc_1 <- cuminc_df$surv_prob[cuminc_df$trt == 1]
    rr <- cuminc_1/cuminc_0
    ve <- 1 - rr

    estimates <- cbind(cuminc_0, cuminc_1,
                       "risk_ratio" = rr,
                       "vaccine_effectiveness" = ve)
    rownames(estimates) <- eval_times


    out <- list(estimates = estimates)

    if(keep_models){
        out$models <- models
    }

    return(out)
}




