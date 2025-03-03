#' Compute Wald or percentile bootstrapped confidence interval
#'
#' @description
#' Internally calls [compute_wald_ci()] and/or [compute_percentile_ci()]
#' depending on the value of `ci_type`
#'
#'
#' @param x A numeric vector of point estimates
#' @param boot_x A matrix of bootstrapped estimates where the rows are the
#' bootstrap iterations and the columns are the time points of interest.
#' @param ci_type Character string indicating which type of confidence interval
#'   to return ("wald", "percentile", "both")
#' @param alpha  Significance level used to compute confidence intervals.
#'   Confidence intervals have nominal level `1 - alpha`.
#' @param transform If `ci_type = "wald"`, a character string indicating the scale on which
#' to compute Wald confidence intervals that are then transformed back to the original scale.
#' Options are `logit` for the transformation log(x/(1-x)) or `log_ve` for the transformation
#' log(1-x).
#' @param z_star If `ci_type = "wald"`, a specific critical value used to
#' to compute Wald confidence intervals (assumed to be positive). If used, `alpha` argument is ignored.
#'
#' @return A matrix containing the lower and upper confidence intervals and
#' relevant bootstrap standard errors.
#' @export
#'
compute_boot_ci <- function(x, boot_x, ci_type, alpha = .05, transform = NULL, z_star = NULL){
    #note boot_x may contains NAs. not able to subset because need to keep matrix structure
    boot_sd <- apply(boot_x, 2, \(x) stats::sd(x, na.rm = TRUE))

    if(ci_type == "wald"){
        ci_out <- compute_wald_ci(x, boot_x, alpha, transform, z_star)

    }else if(ci_type == "percentile"){
        ci_out <- compute_percentile_ci(boot_x, alpha)

    }else if(ci_type == "both"){
        ci_wald <- compute_wald_ci(x, boot_x, alpha, transform, z_star)
        ci_percentile <- compute_percentile_ci(boot_x, alpha)
        ci_out <- cbind(ci_wald, ci_percentile)
    }

    cbind(ci_out, boot_sd)
}


#' @rdname compute_boot_ci
compute_wald_ci <- function(x, boot_x,  alpha = .05, transform, z_star = NULL){
    if(!is.null(z_star)){
        z_crit <- z_star
    }else{
        z_crit <-  stats::qnorm(1 - alpha/2)
    }
    if(transform == "logit"){

        x_logodds  <- logit(x)
        boot_logodds <- logit(boot_x)
        boot_sd <- apply(boot_logodds, 2, \(x) stats::sd(x[is.finite(x)], na.rm = TRUE))
        boot_n <- apply(boot_logodds, 2, \(x) sum(is.finite(x)))
        lower <- stats::plogis(x_logodds - z_crit*boot_sd)
        upper <- stats::plogis(x_logodds + z_crit*boot_sd)
    }else if(transform == "log_ve"){

        x_log <- log_ve(x)
        boot_log <- log_ve(boot_x)
        boot_sd <- apply(boot_log, 2, \(x) stats::sd(x[is.finite(x)], na.rm = TRUE))
        boot_n <- apply(boot_log, 2, \(x) sum(is.finite(x)))
        lower <- -1*(exp(x_log + z_crit*boot_sd) - 1)
        upper <-  -1*(exp(x_log - z_crit*boot_sd) - 1)
    }
    wald_ci <- cbind(lower, upper, boot_sd, wald_n = boot_n)
    colnames(wald_ci) <- c("wald_lower", "wald_upper","wald_sd", "wald_n")
    wald_ci
}

logit <- function(x){log(x/(1-x))}
log_ve <- function(x){log(1-x)}


#' @rdname compute_boot_ci
compute_percentile_ci <- function(boot_x, alpha = .05){
    lower <- apply(boot_x, 2, \(x) stats::quantile(x, alpha/2, na.rm = TRUE))
    upper <- apply(boot_x, 2, \(x) stats::quantile(x, 1 - alpha/2, na.rm = TRUE))
    percentile_ci <- cbind(lower, upper)
    colnames(percentile_ci) <- c("percentile_lower", "percentile_upper")
    percentile_ci
}
