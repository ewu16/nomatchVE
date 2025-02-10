#' Compute a matching-based estimator of VE with confidence intervals
#'
#' @description This function is the main function for computing a matching-based estimator.
#'
#' @inheritParams get_one_matching_ve
#' @inheritParams obsve
# @param matched_data A data frame for the matched cohort
# @param outcome_name Character string representing the original outcome variable in  `matched_data`
# @param event_name Character string representing the original event variable in `matched_data`
# @param trt_name Character string representing the original vaccination status variable in `matched_data`
# @param time_name Character string representing the original vaccination time variable in `matched_data`
# @param adjust_vars  A character vector containing the names of variables in
#   `matched_data` to adjust for.
# @param times  A numeric vector containing the times after vaccination at which to evaluate VE
# @param censor_time Numeric representing the time after vaccination at which individuals should be censored.
# Recommended to set this to `max(times)`
# @param pair_censoring
# @param separate
# @param ci_type
# @param limit_type
#' @param data If limit_type = "limit", the original data used to constructed the matched cohort
#' @param id_name  Character string representing the individual identifier variable  in `data`
#' @param matching_vars A character vector containing the names of variables in `data` to match on.
# @param n_boot
# @param alpha
# @param return_models
# @param return_boot
#'
#' @return A list containing the following:
#' \describe{
#'  \item{estimates}{A list of matrices of the estimates at each timepoint. Rows of
#'  each matrix are the terms "risk_0", "risk_1", "ve". Columns of each matrix
#'  gives the point estimate and confidence intervals at the specified timepoint.}
#'  \item{times}{The timepoints at which VE was evaluated}
#'  \item{n_success_boot}{The number of bootstraps for all the timepoints. TODO: inconsistent with [timepoints()].}
#'  \item{models}{If `return_models = TRUE`, the models used to fit the original point estimates}
#'  \item{boot_samples_list}{If `return_boot = TRUE`, a list of matrices containing the bootstrap estimates for each timepoint}
#' }
#' @export
#'
matching_ve <- function(matched_data,
                        outcome_name,
                        event_name,
                        trt_name,
                        time_name,
                        adjust,
                        times,
                        censor_time,
                        tau,
                        pair_censoring = TRUE,
                        separate = FALSE,
                        ci_type = "wald",
                        limit_type = "fixed",
                        data = NULL,
                        id_name = NULL,
                        matching_vars = NULL,
                        n_boot = 0,
                        alpha = 0.05,
                        return_models = TRUE,
                        return_boot = TRUE){

    # --------------------------------------------------------------------------
    # 1 - Get original estimate
    # --------------------------------------------------------------------------
    original <- get_one_matching_ve(matched_data = matched_data,
                                     outcome_name = outcome_name,
                                     event_name = event_name,
                                     trt_name = trt_name,
                                     time_name = time_name,
                                     adjust = adjust,
                                     times = times,
                                     censor_time = censor_time,
                                     tau = tau,
                                     pair_censoring = pair_censoring,
                                     separate = separate)


    # --------------------------------------------------------------------------
    # 2 - Get bootstrap CI
    # --------------------------------------------------------------------------
    boot_inference <- estimate_matching_ci(matched_data = matched_data,
                                           outcome_name = outcome_name,
                                           event_name = event_name,
                                           trt_name = trt_name,
                                           time_name = time_name,
                                           adjust = adjust,
                                           times = times,
                                           censor_time = censor_time,
                                           tau = tau,
                                           pair_censoring = pair_censoring,
                                           separate = separate,
                                           ci_type = ci_type,
                                           limit_type = limit_type,
                                           data = data,
                                           id_name = id_name,
                                           matching_vars = matching_vars,
                                           n_boot = n_boot,
                                           pt_est = original$estimates,
                                           alpha = alpha,
                                           return_boot = return_boot)

    # --------------------------------------------------------------------------
    # 3 - Format and return results
    # --------------------------------------------------------------------------
    estimates_list <- lapply(seq_along(times), \(i){
        estimate <- c(risk_0 = unname(original$estimates$risk_0[i]),
                      risk_1 = unname(original$estimates$risk_1[i]),
                      ve = unname(original$estimates$ve[i]))
        cbind(estimate,
              boot_inference$ci_estimates_list[[i]])
    })
    names(estimates_list) <- times

    out <- list(estimates = estimates_list,
                times = times,
                n_success_boot = boot_inference$n_success_boot)


    if(return_models){
        out$models <- original$models
    }
    if(return_boot){
        out$boot_samples_list <- boot_inference$boot_samples_list
    }

    return(out)
}


#' Estimate bootstrapped confidence intervals for matching-based estimator
#'
#' @description This function computes Wald and percentile bootstrapped confidence
#' intervals for matching-based VE estimates
#'
#' @inheritParams matching_ve
#' @param pt_est A list of matrices of the estimates for the original data at each timepoint.
#' Rows of each matrix are the terms "risk_0", "risk_1", "ve". Columns of each matrix
#' gives the point estimate and confidence intervals at the specified timepoint.
#'
#' @return A list containing the following:
#' \describe{
#'  \item{ci_estimates_list}{A list of matrices for each timepoint containing the
#'  estimated confidence intervals. Rows of  each matrix are the terms "risk_0", "risk_1", "ve"
#'  Columns of each matrix confidence intervals at the specified timepoint.}
#'  \item{n_success_boot}{The number of bootstraps for all the timepoints. TODO: inconsistent with [timepoints()].}
#'  \item{boot_samples_list}{If `return_boot = TRUE`, a list of matrices containing
#'  the bootstrap estimates for each timepoint.}
#' }
#' @export
#'
estimate_matching_ci <- function(matched_data,
                        outcome_name = outcome_name,
                        event_name,
                        trt_name,
                        time_name,
                        adjust,
                        times,
                        censor_time,
                        tau,
                        pair_censoring = pair_censoring,
                        separate = FALSE,
                        ci_type = "wald",
                        limit_type = "fixed",
                        data = NULL,
                        id_name = NULL,
                        matching_vars = NULL,
                        n_boot = 0,
                        pt_est = NULL,
                        alpha = 0.05,
                        return_boot){
    if(n_boot == 0){
        return(NULL)
    }

    # --------------------------------------------------------------------------
    # 1. Get bootstrap samples
    # --------------------------------------------------------------------------
    cat("Bootstrapping...\n")
    start <- Sys.time()
    empty_mat <- matrix(NA, nrow = n_boot, ncol = length(times),
                        dimnames = list(NULL, times))
    risk_0_mat <- risk_1_mat <- ve_mat <- empty_mat

    for(i in 1:n_boot){
        #set.seed(i)
        if(i %% 50 == 0) print(i)
        boot_estimates <- one_boot_matching(matched_data = matched_data,
                                            outcome_name = outcome_name,
                                            event_name = event_name,
                                            trt_name = trt_name,
                                            time_name = time_name,
                                            adjust = adjust,
                                            times = times,
                                            censor_time = censor_time,
                                            tau = tau,
                                            pair_censoring = pair_censoring,
                                            separate = separate,
                                            limit_type = limit_type,
                                            data = data,
                                            id_name = id_name,
                                            matching_vars = matching_vars,
                                            return_boot = return_boot)


        risk_0_mat[i,] <- boot_estimates$risk_0
        risk_1_mat[i,] <- boot_estimates$risk_1
        ve_mat[i,] <- boot_estimates$ve
    }

    #Boostrap run time
    end <- Sys.time()
    print(end - start)

    # Clean bootstrap samples
    mat_list <- list(risk_0_mat, risk_1_mat, ve_mat)
    check_boot <- lapply(mat_list, \(mat) which(apply(mat, 1, \(x) any(is.na(x)))))
    bad_inds <- unique(unlist(check_boot))
    good_inds <- setdiff(1:n_boot, bad_inds)
    risk_0_mat_clean <- risk_0_mat[good_inds,, drop = FALSE]
    risk_1_mat_clean <- risk_1_mat[good_inds,, drop = FALSE]
    ve_mat_clean <- ve_mat[good_inds,, drop = FALSE]
    n_success_boot <- length(good_inds)
    # print(good_inds)
    # print(risk_0_mat_clean)
    # str(risk_0_mat_clean)

    # --------------------------------------------------------------------------
    # 2. Compute confidence intervals
    # --------------------------------------------------------------------------

    #Temporarily - convert structure to match obsve functions
    boot_samples_list <- lapply(seq_along(times), \(j){
        boot_samples <- matrix(nrow = n_success_boot, ncol = 3,
                               dimnames = list(NULL, c("psi_bar_0", "psi_bar_1", "ve")))
        boot_samples[, 1] <- risk_0_mat_clean[,j]
        boot_samples[, 2] <- risk_1_mat_clean[,j]
        boot_samples[, 3] <- ve_mat_clean[,j]
        boot_samples
    })
    pt_est_list <- lapply(seq_along(times), \(j){
        c(psi_bar_0 = unname(pt_est$risk_0[j]),
          psi_bar_1 = unname(pt_est$risk_1[j]),
          ve = unname(pt_est$ve[j]))
    })

    names(boot_samples_list) <- names(pt_est_list) <- times



    # Compute confidence interval
    if(ci_type == "wald"){
        ci_list <- lapply(seq_along(times), \(i){
            cbind(compute_wald_ci(pt_est_list[[i]], boot_samples_list[[i]], alpha),
                  sd = apply(boot_samples_list[[i]], 2, stats::sd))
            })
    }else if(ci_type == "percentile"){
        ci_list <- lapply(seq_along(times), \(i){
            cbind(compute_percentile_ci(boot_samples_list[[i]], alpha),
                  sd = apply(boot_samples_list[[i]], 2, stats::sd))
        })

    }else if(ci_type == "both"){
        ci_list <- lapply(seq_along(times), \(i){
            cbind(compute_wald_ci(pt_est_list[[i]], boot_samples_list[[i]], alpha),
                  compute_percentile_ci(boot_samples_list[[i]], alpha),
                  boot_sd = apply(boot_samples_list[[i]], 2, stats::sd))
        })
    }
    names(ci_list) <- times

    # --------------------------------------------------------------------------
    # 3. Return confidence intervals and additional bootstrap information
    # --------------------------------------------------------------------------

    out <- list(ci_estimates_list = ci_list,
                n_success_boot = n_success_boot)
    if(return_boot){
        out$boot_samples_list <- boot_samples_list
    }

    return(out)

}


#'  Compute one bootstrap replicate of matching-based VE point estimate
#'
#' @inheritParams matching_ve
#'
#' @return A list with three named elements, "risk_0", "risk_1", "ve". Each
#' element is a vector containing the estimates for the given term at each time point in `times`
#' @export
#'
one_boot_matching <- function(matched_data,
                              outcome_name,
                              event_name,
                              trt_name,
                              time_name,
                              adjust,
                              times,
                              censor_time,
                              tau,
                              pair_censoring = TRUE,
                              separate = FALSE,
                              limit_type = "fixed",
                              data = NULL,
                              id_name = NULL,
                              matching_vars = NULL,
                              return_boot = TRUE){

    # --------------------------------------------------------------------------
    # 1. Create bootstrapped sample(s)
    # --------------------------------------------------------------------------
    if(limit_type == "fixed"){
        #bootstrap from fixed matched cohort
        boot_matched_ids <- sample(unique(matched_data$pair_id), replace = TRUE)
        boot_matched_inds <- as.vector(sapply(boot_matched_ids, \(pair_id) which(matched_data$pair_id == pair_id)))
        boot_matched_data <- matched_data[boot_matched_inds,]
        boot_matched_data$pair_id <- rep(1:length(boot_matched_ids), each = 2)
    }else if(limit_type == "limit"){
        stopifnot("Need to provide original data for limit matching confidence intervals" =
                  !is.null(data))
        boot_inds <- sample(1:nrow(data), replace = TRUE)
        boot_data <-  data[boot_inds,]
        boot_matched_data <- match_rolling_cohort(data = boot_data,
                                                  outcome_name = outcome_name,
                                                  trt_name = trt_name,
                                                  time_name = time_name,
                                                  id_name = id_name,
                                                  matching_vars = matching_vars)[[1]]
    }

    # --------------------------------------------------------------------------
    # 2. Compute VE for bootstrapped data
    # --------------------------------------------------------------------------
    boot_matching_ve <- get_one_matching_ve(matched_data = boot_matched_data,
                                             outcome_name = outcome_name,
                                             event_name = event_name,
                                             trt_name = trt_name,
                                             time_name = time_name,
                                             adjust = adjust,
                                             times = times,
                                             censor_time = censor_time,
                                             tau = tau,
                                             pair_censoring = pair_censoring,
                                             separate = separate,
                                             return_models = FALSE)

    boot_matching_ve$estimates
}




