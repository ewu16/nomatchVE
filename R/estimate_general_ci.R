#' Helper to estimate generic bootstrapped confidence interval
#'
#' @description This function computes Wald and percentile bootstrapped confidence
#' intervals for the proposed VE estimator.
#'
#' @inheritParams obsve
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
#'  \item{boot_samples}{If `return_boot = TRUE`, a list of matrices for each term that contain the bootstrap estimates where the rows are the bootstrap iterations and
#'  the columns are the time points.}
#' }
#' @export
#'
#
estimate_general_ci <- function(one_boot_function_name,
                                one_boot_args,
                                ci_type,
                                n_boot,
                                pt_est = NULL,
                                alpha = 0.05,
                                return_boot = TRUE,
                                n_cores = 1){

    if(n_boot == 0){
        return(NULL)
    }

    # --------------------------------------------------------------------------
    # 1. Run bootstrap
    # --------------------------------------------------------------------------
    cat("Bootstrapping...\n")
    start <- Sys.time()
    times <- one_boot_args$times

    boot_list <- parallel::mclapply(1:n_boot, \(i){
        set.seed(i) #for debugging
        print(i)
        boot_out <- tryCatch({
            output <- capture.output(boot_estimates <- do.call(one_boot_function_name, one_boot_args))
            #using c() here to get a wide vector of estimates to make it easier to create the matrices for each term
            result <- c( boot_estimates[, 1],  boot_estimates[, 2],  boot_estimates[, 3])
            list(result = result ,
                 output = output)
        }, error = function(e){
            list(output = capture.output(print(e)))
        })
    }, mc.cores = n_cores)
    names(boot_list) <- seq_along(boot_list)

    boot_result_list <- lapply(boot_list, \(x) x$result)
    boot_output_list <-  lapply(boot_list, \(x) x$output)

    error_inds <-  which(sapply(boot_result_list, is.null))
    boot_error_list <-  boot_output_list[error_inds]
    success_inds <- setdiff(seq_along(boot_result_list), error_inds)
    boot_success_list <- boot_result_list[success_inds]

    boot_mat <- do.call("rbind", boot_success_list)
    risk_0_mat <- boot_mat[,1:length(times), drop = FALSE]
    risk_1_mat <- boot_mat[,(length(times) + 1): (2*length(times)), drop = FALSE]
    ve_mat <- boot_mat[,(2*length(times) + 1):(3*length(times)), drop = FALSE]
    colnames(risk_0_mat) <- colnames(risk_1_mat) <- colnames(ve_mat) <- times

    #Boostrap run time
    end <- Sys.time()
    print(end - start)

    #if any bootstrap estimate is missing for a specific timepoint, set
    # estimates of other terms in that boostrap to missing
    is_na_mat <- (is.na(risk_0_mat) + is.na(risk_1_mat) + is.na(ve_mat)) > 0
    risk_0_mat[is_na_mat] <- NA
    risk_1_mat[is_na_mat] <- NA
    ve_mat[is_na_mat] <- NA
    n_success_boot <-colSums(!is_na_mat)

    na_inds <- which(rowSums(is_na_mat) > 0)
    boot_na_list <- boot_output_list[na_inds]
    #names(boot_na_list) <- na_inds

    # --------------------------------------------------------------------------
    # 2. Compute confidence intervals
    # --------------------------------------------------------------------------
    risk_0_ci <- compute_boot_ci(x = pt_est[,1],
                                 boot_x = risk_0_mat,
                                 ci_type = ci_type,
                                 alpha = alpha,
                                 transform = "logit")

    risk_1_ci <- compute_boot_ci(x = pt_est[,2],
                                 boot_x = risk_1_mat,
                                 ci_type = ci_type,
                                 alpha = alpha,
                                 transform = "logit")

    ve_ci <- compute_boot_ci(x = pt_est[,3],
                             boot_x = ve_mat,
                             ci_type = ci_type,
                             alpha = alpha,
                             transform = "log_ve")

    # --------------------------------------------------------------------------
    # 3. Return confidence intervals and additional bootstrap information
    # --------------------------------------------------------------------------
    ci_estimates <- list(risk_0 = risk_0_ci,
                         risk_1 = risk_1_ci,
                         ve = ve_ci)

    out <- list(ci_estimates = ci_estimates,
                n_success_boot = n_success_boot,
                boot_error_list = boot_error_list,
                boot_na_list = boot_na_list)

    if(return_boot){
        boot_samples <- list(risk_0 = risk_0_mat,
                             risk_1 = risk_1_mat,
                             ve = ve_mat)
        out$boot_samples <- boot_samples
        out$one_boot_args <- one_boot_args
    }

    return(out)
}

