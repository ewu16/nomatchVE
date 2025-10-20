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
    if(object$effect == "vaccine_effectiveness"){
        title <- "Vaccine Effectiveness Estimates"
        effect_term <- "vaccine_effectiveness"
    }else if(object$effect == "risk_ratio"){
        title <- "Risk Ratio Estimates"
        effect_term <- "risk_ratio"
    }
    cat("\n", title, "\n")
    cat(strrep("=", 50), "\n")

    # Key identifying info
    #cat("Method:", object$method, "\n")
    cat("Call:", paste(deparse(object$call), collapse = "\n"), "\n")

    # The main result
    cat("\nResult:\n")
    effect <- setNames(list(object$estimates[[effect_term]]), effect_term)
    effect_df <- estimates_to_df(effect)

    display_cols <- c("t0", "estimate",
                      grep("lower$", names(effect_df), value = TRUE),
                      grep("upper$", names(effect_df), value = TRUE))

    ci_level <- (1-object$alpha)*100
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

    print(head(effect_df[display_labels], 10), digits = digits)

    # Print how many rows remain
    remaining <- nrow(effect_df) - 10
    if(remaining > 0){
        cat("\nRemaining rows:", remaining, "\n")
    }

    # Hint for more info
    cat("\nUse summary() for more details\n")
    cat("Use plot() to visualize results\n")


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
    cat("Evaluation eval_times:   ", paste(object$eval_times, collapse = ", "), "\n")
    cat("Tau (delay period): ", object$tau, "\n")


    if (length(object$covariates) > 0) {
        cat("Adjusted for:       ", paste(object$covariates, collapse = ", "), "\n")
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


    # cat("N total:")
    # cat("Unexposed:")
    # cat("Exposed:")
    # cat("Exposed remaining at-risk tau days after exposure:")




    # ---- Main Results ----
    cat("\n")
    cat(strrep("-", 70), "\n")
    cat(title, "\n")
    cat(strrep("-", 70), "\n")
    print(head(round(object$estimates[[effect_term]], digits)))

    # Print how many rows remain
    remaining <- nrow(object$estimates[[effect_term]]) - 6
    if(remaining > 0){
        cat("\nRemaining rows:", remaining, "\n")
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
#' @param object An object of class `vefit` created by [nomatchVE()] or [matching_ve()].
#' @param ci_type Character string specifying the type of confidence interval bands to plot.
#'   One of "wald", "percentile", "simul", or "none". Must choose a `ci_type` whose lower
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
plot.vefit <- function(object, ci_type = NULL, color = "#0072B2") {

    ci_type <- validate_confint_type(object, ci_type)
    plot_data <- estimates_to_df(object$estimates)
    plot_data$method <- "dummy"
    alpha <- object$alpha
    plot_ve_panel(plot_data, ci_type, alpha, colors = color,
                  trt_0_label = paste(object$exposure, "= 0"),
                  trt_1_label = paste(object$exposure, "= 1"))
}

#' Check if object is a vefit
#' @keywords internal
#' @noRd
is.vefit <- function(object) {
    inherits(object, "vefit")
}

