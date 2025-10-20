#' Helper to estimate generic bootstrapped confidence interval
#'
#' @description This function computes Wald and percentile bootstrapped confidence
#' intervals for the proposed VE estimator.
#'
#' @inheritParams nomatchVE
#' @param one_boot_function_name Character string naming the function that computes
#' one bootstrap iteration of estimates
#' @param one_boot_args List of arguments to pass to the function that computes
#' one bootstrap iteration of estimations
#' @param  pt_est A matrix of estimates. The columns of the matrix are the cumulative
#' incidence/VE terms and the rows are the requested time points for evaluation.
#'
#' @return A list containing the following:
#' \describe{
#'  \item{ci_estimates}{A list of matrices containing the lower and upper confidence interval bounds and
#'  bootstrap standard error for each term. When `ci_type = "wald"`, the bootstrapped standard errors
#'  on the transformed scale are also included.}
#' \item{n_success_boot}{A numeric vector of the number of successful bootstrap samples for each time point.(Success bootstrap samples are
#'  those that result in non-missing valid point estimates.)}
#'  \item{boot_samples}{If `keep_boot_samples = TRUE`, a list of matrices for each term that contain the bootstrap estimates where the rows are the bootstrap iterations and
#'  the columns are the time points.}
#' }
#' @keywords internal
#'
#
estimate_bootstrap_ci <- function(one_boot_function_name,
                                  one_boot_args,
                                  ci_type,
                                  boot_reps,
                                  pt_est = NULL,
                                  alpha = 0.05,
                                  keep_boot_samples = TRUE,
                                  n_cores = 1){

    if(boot_reps == 0){
        return(NULL)
    }

    # --------------------------------------------------------------------------
    # 1. Run bootstrap
    # --------------------------------------------------------------------------
    cat("Bootstrapping", boot_reps, "samples...\n")
    start <- Sys.time()
    eval_times <- one_boot_args$eval_times

    boot_list <- parallel::mclapply(1:boot_reps, \(i){
        set.seed(i)
        boot_out <- tryCatch({
            output <- utils::capture.output(boot_estimates <- do.call(one_boot_function_name, one_boot_args))
            result <- as.vector(boot_estimates)
            list(result = result ,
                 output = output)
        }, error = function(e){
            list(output = utils::capture.output(print(e)))
        })
    }, mc.cores = n_cores)

    #Bootstrap run time
    end <- Sys.time()
    print(end - start)

    names(boot_list) <- seq_along(boot_list)

    boot_result_list <- lapply(boot_list, \(x) x$result)
    boot_output_list <-  lapply(boot_list, \(x) x$output)

    error_inds <-  which(sapply(boot_result_list, is.null))
    boot_error_list <-  boot_output_list[error_inds]
    success_inds <- setdiff(seq_along(boot_result_list), error_inds)
    boot_success_list <- boot_result_list[success_inds]

    if(length(boot_success_list) == 0){
        stop("All bootstrap replicates encountered an error, e.g.\n",
             boot_error_list[[1]])

    }

    boot_mat <- do.call("rbind", boot_success_list)
    cuminc_0_mat <- boot_mat[,1:length(eval_times), drop = FALSE]
    cuminc_1_mat <- boot_mat[,(length(eval_times) + 1): (2*length(eval_times)), drop = FALSE]
    rr_mat <- boot_mat[,(2*length(eval_times) + 1):(3*length(eval_times)), drop = FALSE]
    ve_mat <- boot_mat[,(3*length(eval_times) + 1):(4*length(eval_times)), drop = FALSE]
    colnames(cuminc_0_mat) <- colnames(cuminc_1_mat) <- colnames(rr_mat) <- colnames(ve_mat) <- eval_times



    #if any bootstrap estimate is missing for a specific timepoint, set
    # estimates of other terms in that boostrap to missing
    is_na_mat <- (is.na(cuminc_0_mat) + is.na(cuminc_1_mat) + is.na(rr_mat) + is.na(ve_mat)) > 0
    cuminc_0_mat[is_na_mat] <- NA
    cuminc_1_mat[is_na_mat] <- NA
    rr_mat[is_na_mat] <- NA
    ve_mat[is_na_mat] <- NA
    n_success_boot <-colSums(!is_na_mat)

    na_inds <- which(rowSums(is_na_mat) > 0)
    boot_na_list <- boot_output_list[na_inds]
    #names(boot_na_list) <- na_inds

    # --------------------------------------------------------------------------
    # 2. Compute confidence intervals
    # --------------------------------------------------------------------------
    cuminc_0_ci <- compute_boot_ci(x = pt_est[, "cuminc_0"],
                                 boot_x = cuminc_0_mat,
                                 ci_type = ci_type,
                                 alpha = alpha,
                                 transform = "logit")

    cuminc_1_ci <- compute_boot_ci(x = pt_est[, "cuminc_1"],
                                 boot_x = cuminc_1_mat,
                                 ci_type = ci_type,
                                 alpha = alpha,
                                 transform = "logit")

    rr_ci <- compute_boot_ci(x = pt_est[, "risk_ratio"],
                                   boot_x = rr_mat,
                                   ci_type = ci_type,
                                   alpha = alpha,
                                   transform = "log_rr")

    ve_ci <- compute_boot_ci(x = pt_est[, "vaccine_effectiveness"],
                             boot_x = ve_mat,
                             ci_type = ci_type,
                             alpha = alpha,
                             transform = "log_ve")

    # --------------------------------------------------------------------------
    # 3. Return confidence intervals and additional bootstrap information
    # --------------------------------------------------------------------------
    ci_estimates <- list(cuminc_0 = cuminc_0_ci,
                         cuminc_1 = cuminc_1_ci,
                         risk_ratio = rr_ci,
                         vaccine_effectiveness = ve_ci)

    out <- list(ci_estimates = ci_estimates,
                n_success_boot = n_success_boot,
                boot_error_list = boot_error_list,
                boot_na_list = boot_na_list)

    if(keep_boot_samples){
        boot_samples <- list(cuminc_0 = cuminc_0_mat,
                             cuminc_1 = cuminc_1_mat,
                             risk_ratio = rr_mat,
                             vaccine_effectiveness = ve_mat)
        out$boot_samples <- boot_samples
        out$one_boot_args <- one_boot_args
    }

    return(out)
}

