#' Compare two vefit objects
#'
#' @description
#' Plot cumulative incidence and VE estimates for two
#' different methods using colors to distinguish methods
#'
#' @param fit1 A vefit object (typically from \code{\link{matching_ve}})
#' @param fit2 A vefit object (typically from \code{\link{nomatchVE}})
#' @param labels Character vector of length 2 providing labels for the two methods.
#'  Default is \code{c("Method 1", "Method 2")}.
#' @param confint_type Character string specifying the type of confidence interval to plot.
#'   One of "wald", "percentile", or "simul". If \code{NULL} (default), uses the
#'   CI type from \code{fit1}. If the object has \code{confint_type = "both"},
#'   defaults to "wald".
#' @param colors Character vector of length 2 providing colors for the two methods.
#'   Default is \code{c("#F8766D", "#00BFC4")} (ggplot2's default red and cyan).
#'
#' @return A ggplot2 object with three faceted panels showing cumulative incidence
#'   and VE estimates for both methods.
#'
#' @details
#' Both \code{fit1} and \code{fit2} must have the same significance level (alpha).
#' The function will stop with an error if alphas differ.
#'
#' For cumulative incidence panels, y-axis limits are shared across methods to
#' facilitate comparison. The VE panel uses free y-axis scaling.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit both methods
#' fit_nomatch <- nomatchVE(data = simdata, ...)
#' fit_match <- matching_ve(matched_data = matched_data, ...)
#' # Compare with custom labels
#' compare_ve_fits(
#'   fit_match,
#'   fit_nomatch,
#'   labels = c("Matching", "G-computation"))
#' }
compare_ve_fits <- function(fit1,
                            fit2,
                            labels = c("Method 1", "Method 2"),
                            confint_type = NULL,
                            colors = c("#F8766D", "#00BFC4")) {

    # Validation
    if (!is.vefit(fit1) || !is.vefit(fit2)) {
        stop("Both fit1 and fit2 must be vefit objects", call. = FALSE)
    }

    if (length(labels) != 2 || length(colors) != 2) {
        stop("labels and colors must be vectors of length 2", call. = FALSE)
    }

    confint_type <- validate_confint_type(fit1, confint_type)
    confint_type <- validate_confint_type(fit2, confint_type)

    alpha1 <- fit1$alpha
    alpha2 <- fit2$alpha
    if(alpha1 != alpha2){
        stop("Both fit1 and fit2 must have same alpha level")
    }

    # Prepare data from both fits
    data1 <- estimates_to_df(fit1$estimates)
    data2 <- estimates_to_df(fit2$estimates)

    # Add method labels and rename CI columns
    data1$method <- labels[1]
    data2$method <- labels[2]

    # Combine
    plot_data <- rbind(data1, data2)
    plot_data$method <- factor(plot_data$method, levels = labels)

    # Call internal plotting function
    plot_ve_panel(
        plot_data = plot_data,
        confint_type = confint_type,
        alpha = alpha1,
        colors = colors
    )
}

