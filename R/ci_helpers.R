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
#' @return A matrix containing the lower and upper confidence intervals
#'
#' @keywords internal
#' @noRd
#'
compute_boot_ci <- function(x, boot_x, ci_type, alpha = .05, transform = NULL, z_star = NULL){
    if(ci_type == "wald"){
        ci_out <- compute_wald_ci(x, boot_x, alpha, transform, z_star)

    }else if(ci_type == "percentile"){
        ci_out <- compute_percentile_ci(boot_x, alpha)

    }else if(ci_type == "both"){
        ci_wald <- compute_wald_ci(x, boot_x, alpha, transform, z_star)
        ci_percentile <- compute_percentile_ci(boot_x, alpha)
        ci_out <- cbind(ci_wald, ci_percentile)
    }

    ci_out
}


#' @rdname compute_boot_ci
#' @keywords internal
#' @noRd
compute_wald_ci <- function(x, boot_x,  alpha = .05, transform, z_star = NULL){

    z <- if (!is.null(z_star)) z_star else stats::qnorm(1 - alpha/2)

    #Define the functions for forward and back transformation
    tf <- switch(transform,
        "logit" = list(
            fwd   = stats::qlogis,
            lower = function(eta, sd) stats::plogis(eta - z * sd),
            upper = function(eta, sd) stats::plogis(eta + z * sd)
        ),
        "log_ve" = list(
            fwd   = function(y) log(1-y),
            lower = function(eta, sd) 1 - exp(eta + z * sd),
            upper = function(eta, sd) 1 - exp(eta - z * sd)
        ),
        "log_rr" = list(
            fwd   = log,
            lower = function(eta, sd) exp(eta - z * sd),
            upper = function(eta, sd) exp(eta + z * sd)
        )
    )

    # transform
    eta_x    <- tf$fwd(x)
    eta_boot <- tf$fwd(boot_x)

    # bootstrap SDs/counts ignoring non-finite draws
    col_sd  <- apply(eta_boot, 2, function(col) stats::sd(col[is.finite(col)], na.rm = TRUE))
    boot_n  <- apply(eta_boot, 2, function(col) sum(is.finite(col)))

    # back-transform CI limits
    lower <- tf$lower(eta_x, col_sd)
    upper <- tf$upper(eta_x, col_sd)

    cbind(wald_lower = lower, wald_upper = upper, wald_n = boot_n)

}


#' @rdname compute_boot_ci
#' @keywords internal
#' @noRd
compute_percentile_ci <- function(boot_x, alpha = .05){
    lower <- apply(boot_x, 2, \(x) stats::quantile(x, alpha/2, na.rm = TRUE))
    upper <- apply(boot_x, 2, \(x) stats::quantile(x, 1 - alpha/2, na.rm = TRUE))
    boot_n <- apply(boot_x, 2, \(x) sum(!is.na(x)))

    cbind(percentile_lower = lower, percentile_upper = upper, percentile_n = boot_n)
}
