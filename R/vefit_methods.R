#' Print method for vaccine effectiveness fits
#'
#' @description
#' Prints a concise summary of vaccine effectiveness estimates from a fitted model.
#'
#' @param object An object of class `vefit` created by [nomatchVE()] or [matching_ve()].
#' @param digits Integer indicating the number of decimal places to display. Default is 3.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' Invisibly returns the input object `object`. Called primarily for its side effect
#' of printing a summary to the console.
#'
#'
#'
#' @export
print.vefit <- function(object, digits = 3, ...) {
    cat("\nVaccine Effectiveness Estimates\n")
    cat(strrep("=", 50), "\n")

    # Key identifying info
    cat("Method:", object$method, "\n")
    cat("Times:", paste(object$times, collapse = ", "), "\n")
    cat("Tau:", object$tau, "\n")

    # The main result - VE estimates
    cat("\nVE Estimates:\n")
    ve_mat <- object$estimates$ve

    # Show estimate and CI if available
    display_cols <- c("estimate",
                      grep("lower$", colnames(ve_mat), value = TRUE),
                      grep("upper$", colnames(ve_mat), value = TRUE))
    display_cols <- display_cols[!is.na(display_cols)]
    print(round(ve_mat[, display_cols, drop = FALSE], digits))


    # Hint for more info
    cat("\n")
    cat("Call summary() for more details\n")
    cat("Call plot() to visualize results\n")


    invisible(object)
}

#' Summary method for vaccine effectiveness fits
#' @description
#' Summarizes how vaccine effectiveness estimates were obtained and displays
#' cumulative incidence and VE across all evaluation time points.
#'
#' @param object An object of class `vefit` created by [nomatchVE()] or [matching_ve()].
#' @param digits Integer indicating the number of decimal places to display. Default is 4.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' Invisibly returns the input object `object`. Called primarily for its side effect
#' of printing a detailed summary to the console.
#'
#' @export
summary.vefit <- function(object, digits = 4, ...) {
    cat("\n")
    cat(strrep("=", 70), "\n")
    cat("Vaccine Effectiveness Analysis Summary\n")
    cat(strrep("=", 70), "\n\n")

    # ---- Key Parameters ----
    cat("Method:             ", object$method, "\n")
    cat("Evaluation times:   ", paste(object$times, collapse = ", "), "\n")
    cat("Tau (delay period): ", object$tau, "\n")
    cat("Censoring time:     ", object$censor_time, "\n")

    if (length(object$adjust_vars) > 0) {
        cat("Adjusted for:       ", paste(object$adjust_vars, collapse = ", "), "\n")
    }

    # ---- Bootstrap Info (if applicable) ----
    if (object$n_boot > 0) {
        cat("\nBootstrap:          ", object$n_boot, "replicates\n")
        cat("Confidence level:   ", (1 - object$alpha) * 100, "%\n")
        cat("Successful samples: ",
            paste(range(object$n_success_boot), collapse = "-"),
            " (range across timepoints)\n")
    } else {
        cat("\nNo bootstrap performed (n_boot = 0)\n")
    }

    # ---- Main Results ----
    cat("\n")
    cat(strrep("-", 70), "\n")
    cat("Cumulative Incidence:", object$trt_name, "= 0\n")
    cat(strrep("-", 70), "\n")
    print(round(object$estimates$cuminc_0, digits))

    cat("\n")
    cat(strrep("-", 70), "\n")
    cat("Cumulative Incidence:", object$trt_name, "= 1\n")
    cat(strrep("-", 70), "\n")
    print(round(object$estimates$cuminc_1, digits))

    cat("\n")
    cat(strrep("-", 70), "\n")
    cat("Vaccine Effectiveness\n")
    cat(strrep("-", 70), "\n")
    print(round(object$estimates$ve, digits))


    cat("\n")
    cat(strrep("=", 70), "\n\n")

    invisible(object)
}

#' Plot method for vefit objects
#'
#' @description
#' Plot cumulative incidence and VE estimates across all evaluation time points.
#'
#'
#' @param object An object of class `vefit` created by [nomatchVE()] or [matching_ve()].
#' @param confint_type Character string specifying the type of confidence interval bands to plot.
#'   One of "wald", "percentile", "simul", or "none". Must choose a `confint_type` whose lower
#'   and upper bounds are already computed in `object$estimates` component of `object`.
#' @param color Aesthetic value to map data values to. Default: `"#0072B2"` (blue)
#'
#' @return a ggpplot2 object with line plot for each term (cumulative incidence, VE)
#'
#' @details
#' For cumulative incidence panels, y-axis limits are shared across methods to
#' facilitate comparison. The VE panel uses free y-axis scaling.
#'
#' @export
plot.vefit <- function(object, confint_type = NULL, color = "#0072B2") {

    confint_type <- validate_confint_type(object, confint_type)
    plot_data <- estimates_to_df(object$estimates)
    plot_data$method <- "dummy"
    alpha <- object$alpha
    plot_ve_panel(plot_data, confint_type, alpha, colors = color)
}

#' Check if object is a vefit
#' @keywords internal
#' @noRd
is.vefit <- function(object) {
    inherits(object, "vefit")
}

