#' Helper function to estimate marginal risk in matched data set
#'
#' @description
#' First, creates the analysis matched data set. Then calls [estimate_marginal_risk()],
#' providing the appropriate outcomes based on the analysis data set.
#'
#' @inheritParams clean_matched_data
#' @inheritParams estimate_marginal_risk
#'
#' @return The object returned by [estimate_marginal_risk()]
#' @export
#'
get_one_matching_ve <- function(matched_data,
                                 outcome_name,
                                 event_name,
                                 trt_name,
                                 time_name,
                                 adjust,
                                 times,
                                 censor_time,
                                 tau,
                                 pair_censoring = TRUE,
                                 separate = TRUE,
                                 return_models = TRUE){

    matched_adata <- clean_matched_data(matched_data = matched_data,
                                        outcome_name = outcome_name,
                                        event_name = event_name,
                                        trt_name = trt_name,
                                        time_name = time_name,
                                        tau = tau,
                                        pair_censoring = pair_censoring)

    estimate_marginal_risk(adata = matched_adata,
                           adata_outcome_name = "T_d",
                           adata_event_name = paste0(event_name, "_d"),
                           adata_trt_name = paste0(trt_name, "_d"),
                           adjust = adjust,
                           times = times,
                           censor_time = censor_time,
                           tau = tau,
                           separate = separate,
                           return_models = return_models)
}





#' Get the marginal cumulative incidence in the treated and untreated
#' groups based on a Cox model.
#'
#' @description This function uses a Cox model to predict the marginal cumulative
#' incidence of an endpoint `t0` days after vaccination. The marginal
#' cumulative incidence is computed by predicting individual survival probabilities
#' and averaging these predictions over everyone in the treatment group.
#'
#' @param adata A data fram that represents the analysis data set of a clinical trial.
#' @param adata_outcome_name Character string specifying the time to event variable in `adata`. The
#' time should be the time to event from vaccination/matched index date
#' @param adata_event_name Character string specifying the event variable in `adata`
#' @param adata_trt_name Character string specifying the treatment variable in `adata`
#' @param adjust A character vector  of variables in `adata` to adjust for, or a formula object.
#' Do not include treatment variable.
#' @param tau waiting period
#' @param times A numeric vector of the time points at which to return cumulative incidence estimates
#' @param censor_time Time at which observations should be censored. Typically the `max(times)`.
#' @param separate If `TRUE`, models are fit in the treated and untreated groups separately
#' @param return_models If `TRUE`, fitted survival mdoels are returned.

#'
#' @return A list containing the following:
#' \describe{
#' \item{estimates}{A list with three named elements, "risk_0", "risk_1", "ve". Each
#' element is a vector containing the estimates for the given term at each time point in `times`}
#' \item{models}{If `return_models = TRUE`, the models used to compute risk and VE}
#' }
#' @export
#'
estimate_marginal_risk <- function(adata,
                                   adata_outcome_name,
                                   adata_event_name,
                                   adata_trt_name,
                                   adjust,
                                   times,
                                   censor_time,
                                   tau,
                                   separate = FALSE,
                                   return_models = TRUE){
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
    models <- lapply(data_list, \(df){survival::coxph(model_formula, data = df, model = TRUE)})

    # --------------------------------------------------------------------------
    # 2. Make predictions
    # --------------------------------------------------------------------------

    # Do prediction and compute marginal risk
    predict_at_t <- function(model_data, model, t){
        pred_data <- model_data
        pred_data$tt <- t
        surv_prob <- stats::predict(model, newdata = pred_data, type = "survival")
        data.frame(t0 = t, trt = model_data$trt, surv_prob)
    }

    pred_at_times <- function(model_data, model, times){
        pred_list <- lapply(times, \(t) predict_at_t(model_data, model, t))
        names(pred_list) <- times
        do.call("rbind", pred_list)
        rownames(pred_list) <- NULL
        pred_list
    }


    combine_list <- lapply(seq_along(models), \(i){
       pred_at_times(model_data = data_list[[i]],
                     model = models[[i]],
                     times)
    })

    # --------------------------------------------------------------------------
    # 3. Marginalize
    # --------------------------------------------------------------------------


    full_pred <- do.call("rbind", unlist(combine_list, recursive = FALSE))
    rownames(full_pred) <- NULL
    #by using aggregate can keep things simpler for separate and together case
    risk_df <- stats::aggregate(surv_prob ~trt + t0, data = full_pred, FUN =  \(x){1 - mean(x)})


    # --------------------------------------------------------------------------
    # 4. Format result
    # --------------------------------------------------------------------------

    risk_0 <- risk_df$surv_prob[risk_df$trt == 0]
    risk_1 <- risk_df$surv_prob[risk_df$trt == 1]
    ve <- 1 - risk_1/risk_0
    names(ve) <- names(risk_1) <- names(risk_0) <-times

    estimates <- list(risk_0 = risk_0,
                      risk_1 = risk_1,
                      ve = ve)

    out <- list(estimates = estimates)

    if(return_models){
       out$models <- models
    }

    return(out)

}




