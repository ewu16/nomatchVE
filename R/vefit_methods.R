#' Print method for vaccine effectiveness fits
#'
#' @description
#' Prints a concise summary of vaccine effectiveness estimates from a fitted model.
#'
#' @param x An x of class `vefit` created by [nomatchVE()] or [matching_ve()].
#' @param digits Integer indicating the number of decimal places to display. Default is 3.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' Invisibly returns the input x `x`. Called primarily for its side effect
#' of printing a summary to the console.
#'
#'
#'
#' @export
print.vefit <- function(x, digits = 3, ...) {
    if(x$effect == "vaccine_effectiveness"){
        title <- "Vaccine Effectiveness Estimates"
        effect_term <- "vaccine_effectiveness"
    }else if(x$effect == "risk_ratio"){
        title <- "Risk Ratio Estimates"
        effect_term <- "risk_ratio"
    }
    cat("\n", title, "\n")
    cat(strrep("=", 50), "\n")

    # Key identifying info
    #cat("Method:", x$method, "\n")
    cat("Call:", paste(deparse(x$call), collapse = "\n"), "\n")

    # The main result
    cat("\nResult:\n")
    effect <- stats::setNames(list(x$estimates[[effect_term]]), effect_term)
    effect_df <- estimates_to_df(effect)

    display_cols <- c("t0", "estimate",
                      grep("lower$", names(effect_df), value = TRUE),
                      grep("upper$", names(effect_df), value = TRUE))

    ci_level <- (1-x$alpha)*100
    name_map <- c(
        t0                = "Timepoint",
        estimate          = "Estimate",
        wald_lower        = paste0(ci_level, "% Wald CI: Lower"),
        wald_upper        = paste0(ci_level, "% Wald CI: Upper"),
        percentile_lower  = paste0(ci_level, "% Percentile CI: Lower"),
        percentile_upper  = paste0(ci_level, "% Percentile CI: Upper")
    )

    display_labels <- name_map[display_cols]
    names(effect_df)[names(effect_df) %in% display_cols] <- display_labels

    print(utils::head(effect_df[display_labels], 10), digits = digits)

    # Print how many rows remain
    remaining <- nrow(effect_df) - 10
    if(remaining > 0){
        cat("\nRemaining rows:", remaining, "\n")
    }

    # Hint for more info
    cat("\nUse summary() for more details\n")
    cat("Use plot() to visualize results\n")


    invisible(x)
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
    if(object$effect == "vaccine_effectiveness"){
        title <- "Vaccine Effectiveness"
        effect_term <- "vaccine_effectiveness"
    }else if(object$effect == "risk_ratio"){
        title <- "Risk Ratio"
        effect_term <- "risk_ratio"
    }

    cat("\n")
    cat(strrep("=", 70), "\n")
    cat(title, "Analysis Summary\n")
    cat(strrep("=", 70), "\n\n")

    # ---- Key Parameters ----
    cat("Method:             ", object$method, "\n")
    cat("Evaluation times:   ", paste(utils::head(object$eval_times), collapse = ", "),
        ifelse(length(object$eval_times) > 6, ", ...", ""), "\n")
    cat("Tau (delay period): ", object$tau, "\n")

    if(object$method ==  "nomatchVE (G-computation)"){
        if (length(object$covariates) > 0) {
            cat("Adjusted for:       ", paste(object$covariates, collapse = ", "), "\n")
        }
    }



    # ---- Bootstrap Info (if applicable) ----
    if (object$boot_reps > 0) {
        cat("\nBootstrap:          ", object$boot_reps, "replicates\n")
        cat("Confidence level:   ", (1 - object$alpha) * 100, "%\n")
        cat("Successful samples: ",
            paste(range(object$n_success_boot), collapse = "-"),
            " (range across timepoints)\n")
    } else {
        cat("\nNo bootstrap performed (boot_reps = 0)\n")
    }


    if(object$method ==  "nomatchVE (G-computation)"){
        cat("\n")
        cat(strrep("-", 70), "\n")
        cat("Sample:\n")
        cat(strrep("-", 70), "\n")

        descrip <- object$descrip
        cat("N total:", descrip$n, "\n")
        cat("Number of events:", descrip$n_events, "\n")
        cat("\n")
        cat("N exposed:", descrip$n_exposed, "\n")
        cat("N exposed at-risk <tau> days after exposure:", descrip$n_exposed_at_tau, "\n")

        cat("\n")
        cat("Distribution of exposure times among at-risk <tau> days after exposure:\n")
        cat(" Range: ", paste(range(descrip$exposure_times_at_tau), collapse = " - "), "| ")
        cat(" Median (IQR): ",
            stats::median(descrip$exposure_times_at_tau),
            paste0("(", stats::quantile(descrip$exposure_times_at_tau, 0.25), " - ",
                   stats::quantile(descrip$exposure_times_at_tau, 0.75), ")"), "| ")
        cat(" Mean: ",  round(mean(descrip$exposure_times_at_tau), 1))

        # ---- Main Results ----
        cat("\n\n")
        cat(strrep("-", 70), "\n")
        cat("Model for unexposed:\n")
        cat(strrep("-", 70), "\n")
        cat("N =", object$model_0$n, "| Number of events =", object$model_0$nevent, "\n\n")
        print(round(stats::coef(summary(object$model_0)), 3))

        cat("\n")
        cat(strrep("-", 70), "\n")
        cat("Model for exposed:\n")
        cat(strrep("-", 70), "\n")
        cat("N =", object$model_1$n, "| Number of events =", object$model_1$nevent, "\n\n")
        print(round(stats::coef(summary(object$model_1)), 3))

    }else{
        cat("\n")
        cat(strrep("-", 70), "\n")
        cat("Sample:\n")
        cat(strrep("-", 70), "\n")

        descrip <- object$descrip
        cat("N matched:", descrip$n_matched, "\n")
        cat("N matched analysis:",  descrip$n_matched_analysis, "\n",
            "(excludes pairs in which either individual is not at risk\n <tau> days after matching index day)\n\n")
        cat("Number of events in matched analysis:", descrip$n_events, "\n")
        cat("\n")
        cat("Distribution of exposure times in matched analysis:\n")
        cat(" Range: ", paste(range(descrip$exposure_times), collapse = " - "), "| ")
        cat(" Median (IQR): ",
            stats::median(descrip$exposure_times),
            paste0("(", stats::quantile(descrip$exposure_times, 0.25), " - ",
                        stats::quantile(descrip$exposure_times, 0.75), ")"), "| ")
        cat(" Mean: ",  round(mean(descrip$exposure_times), 1))

        cat("\n\n")
        cat(strrep("-", 70), "\n")
        cat("Kaplan Meier for matched analysis:\n")
        cat(strrep("-", 70), "\n\n")
        print(object$models)

    }


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
#' @param x An x of class `vefit` created by [nomatchVE()] or [matching_ve()].
#' @param ci_type Character string specifying the type of confidence interval bands to plot.
#'   One of `"wald", "percentile", "simul"`, or `"none"`. Must choose a `ci_type` whose lower
#'   and upper bounds are already computed in `x$estimates` component of `x`.
#' @param color Aesthetic value to map data values to. Default: `"#0072B2"` (blue)
#' @param ... Additional arguments (currently ignored).
#' @return a ggplot2 x with line plot for each term (cumulative incidence, VE)
#'
#' @details
#' For cumulative incidence panels, y-axis limits are shared across methods to
#' facilitate comparison. The VE panel uses free y-axis scaling.
#'
#' @importFrom rlang .data
#' @export
plot.vefit <- function(x, ci_type = NULL, color = "#0072B2", ...) {

    ci_type <- validate_confint_type(x, ci_type)
    plot_data <- estimates_to_df(x$estimates)
    plot_data$method <- "dummy"
    alpha <- x$alpha
    plot_ve_panel(plot_data, ci_type, alpha, colors = color,
                  trt_0_label = paste(x$exposure, "= 0"),
                  trt_1_label = paste(x$exposure, "= 1"))
}

#' Check if x is a vefit
#' @keywords internal
#' @noRd
is.vefit <- function(x) {
    inherits(x, "vefit")
}

