#' Helper for constructing bootstrapped confidence intervals
#'
#' @description This function computes Wald and percentile bootstrapped
#'   confidence intervals.
#'
#' @inheritParams nomatchVE
#' @param one_boot_function Function that computes one bootstrap iteration and
#'   returns the bootstrap estimates
#' @param one_boot_args List of arguments to pass to `one_boot_function`
#'
#' @param pt_est A matrix of point estimates, with columns corresponding to
#'   cumulative incidence and effect measures, and rows representing evaluation
#'   timepoints. This argument must be provided when `ci_type = "wald"`. Ignored
#'   when `ci_type = "percentile"`.
#'
#' @return A list containing the following:
#' \describe{
#'  \item{ci_estimates}{A list of matrices, one for each term, containing the
#'   lower and upper confidence interval bounds}
#' \item{n_success_boot}{A numeric vector of the number of successful bootstrap
#' samples for each time point.(Success bootstrap samples are
#'  those that result in non-missing valid point estimates.)}
#'  \item{boot_samples}{If `keep_boot_samples = TRUE`, a list of matrices, one
#'  for each term, containing the bootstrap estimates. Rows are the
#'  bootstrap iterations and columns are the time points.}
#' }
#'
#'  When `boot_reps = 0`, returns `NULL`.
#'
#'
#' @keywords internal
#'
estimate_bootstrap_ci <- function(one_boot_function,
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
    if(is.null(pt_est) & ci_type %in% c("wald", "both")){
        stop("Must provide point estimate to construct Wald confidence intervals")
    }

    # --------------------------------------------------------------------------
    # 1. Run bootstrap
    # --------------------------------------------------------------------------
    cat("Bootstrapping", boot_reps, "samples...\n")
    start <- Sys.time()

    boot_list <- parallel::mclapply(1:boot_reps, \(i){
        set.seed(i)
        old <- options(warn = 1); on.exit(options(old), add = TRUE)
        boot_out <- tryCatch({
            output <- utils::capture.output(
                boot_est<- do.call(one_boot_function, one_boot_args),
                type = "message"
            )
            list(result = boot_est,
                 output = output)
        }, error = function(e){
            list(result = NULL,
                 output = utils::capture.output(print(e)))
        })
    }, mc.cores = n_cores)

    #store the bootstrap index number for checking logs
    names(boot_list) <- seq_along(boot_list)

    end <- Sys.time()
    print(end - start)

    # Check for errors
    results  <- lapply(boot_list, \(x) x$result)
    logs     <- lapply(boot_list, \(x) x$output)
    bad      <- which(sapply(results, is.null))
    good     <- setdiff(seq_along(results), bad)
    errors   <- logs[bad]

    if(length(good) == 0){
        stop("All bootstrap replicates encountered an error, e.g.\n ",
             paste(logs[[1]], collapse = "\n"))
    }
    if(length(bad) > 0){
        warning("Some bootstrap replicates encounted an error.\n ",
                "Please check `$n_success_boot` and `$boot_error` for more information.")
    }

    # Compile results
    boot_mat <- do.call(rbind, results[good])
    eval_times <- one_boot_args$eval_times

    boot_samples <- stats::setNames(
        lapply(1:ncol(boot_mat), \(i){
            m <- matrix(as.numeric(boot_mat[,i]), ncol = length(eval_times), byrow = TRUE)
            colnames(m) <- eval_times
            rownames(m) <- 1:nrow(m)
            return(m)
        }),
        colnames(boot_mat)
    )

    # Clean up results
    # if any estimate is missing, set corresponding estimates for all terms to missing
    is_na_mat <- Reduce("+", lapply(boot_samples, is.na)) > 0
    if(any(is_na_mat)){
        warning("Some bootstrap replicates resulted in estimates of NA.\n ",
                "Please check `$n_success_boot` and `$boot_nas` for more information.")
    }
    n_success_boot <-colSums(!is_na_mat)

    boot_samples_clean <- lapply(boot_samples, \(x){
        x[is_na_mat] <- NA
        return(x)
        })

    nas <-logs[which(rowSums(is_na_mat) > 0)]

    # --------------------------------------------------------------------------
    # 2. Compute confidence intervals
    # --------------------------------------------------------------------------
    tr_map <- c(cuminc_0 = "logit",
                cuminc_1 = "logit",
                risk_ratio = "log_rr",
                vaccine_effectiveness = "log_ve")


    ci_estimates <- lapply(names(boot_samples_clean), \(term){
        compute_boot_ci(x = pt_est[, term],
                        boot_x = boot_samples_clean[[term]],
                        ci_type = ci_type,
                        alpha = alpha,
                        transform = tr_map[term])
    })
    names(ci_estimates) <- names(boot_samples_clean)

    # --------------------------------------------------------------------------
    # 3. Return confidence intervals and additional bootstrap information
    # --------------------------------------------------------------------------

    out <- list(ci_estimates = ci_estimates,
                n_success_boot = n_success_boot,
                boot_errors = errors ,
                boot_nas = nas,
                boot_samples = if(keep_boot_samples) boot_samples_clean else NULL
                )

    return(out)
}

