#' Evaluate VE at various time points of interest
#'
#' @param object An object returned from `obsve()`
#' @param times A numeric vector containing the times after vaccination at which to evaluate VE
#'
#' @return The original `obsve` object with the following changes/additions:
#' \describe{
#' \item{estimates}{A list of dataframe/matrix of the estimates at each timepoint}
#' \item{n_success_boot}{A vector of the number of bootstraps for each timepoint}
#' \item{boot_samples}{A list of matrices containing the bootstrap estimates for each timepoint}
#' \item{timepoints}{The timepoints at which VE was evaluated}
#' }
#'
#' @export
#'
timepoints <- function(object, times){

    #Set arguments to call for obsve
    call_list <- as.list(object$call)[-1]

    args_manual <- list(marginalizing_dist = object$gp_list,
                        formula_0 = object$model_0,
                        formula_1 = object$model_1)

    args_keep <-object[which(names(object) %in% names(call_list))]

    args_dummy <- lapply(call_list[which(!names(call_list) %in% c(names(object), names(args_manual)))],
                         \(x) NULL)
    args_list <- c(args_manual, args_keep, args_dummy)

    ##clean up timepoints; also include original timepoint in results
    times_clean <- sort(unique(c(times, object$t0)))
    times_clean <- times_clean[times_clean <= object$t0]
    stopifnot("must choose times <= t0" = length(times_clean) > 0)

    ##create list of estimates
    # estimates_mat <- matrix(NA,
    #                         nrow = length(times_clean),
    #                         ncol = ncol(object$estimates),
    #                         dimnames = list(NULL,names(object$estimates)))
    estimates_list <- vector(mode = "list", length(times_clean))
    n_success_boot <- rep(NA, length(times_clean))
    boot_samples <- vector(mode = "list", length(times_clean))
    names(estimates_list) <-  names(n_success_boot) <- names(boot_samples) <-times_clean

    ####get estimates for each timepoint
    counter <- 0
    for(t in times_clean){
        counter <- counter + 1
        if(t == object$t0){
            my_object <- object

        }else{
            cat("Timepoint:", t, "\n")
            args_list$t0 <- t
            my_object <- do.call("obsve", args = args_list)
        }
        estimates_list[[counter]] <-my_object$estimates
        n_success_boot[counter] <- my_object$n_success_boot
        boot_samples[[counter]] <- my_object$boot_samples
    }


    out  <- object
    out$estimates <- estimates_list
    out$times <- times_clean
    out$n_success_boot <- n_success_boot
    out$boot_samples <- boot_samples

    return(out)
}

