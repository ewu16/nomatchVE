#' Internal plotting function for VE panel plots
#'
#' @description
#' Creates a three-panel plot showing cumulative incidence (no vaccine),
#' cumulative incidence (vaccine), and vaccine effectiveness over time.
#' This function handles both single method and comparison plots.
#'
#' @param plot_data A data frame containing estimates to plot. Must include
#'   columns: `t0` (time), `estimate`, `term` (one of "cuminc_0", "cuminc_1", "ve"),
#'   `<confint_type>_lower`, `<confint_type>_upper`, and `method`
#' @param confint_type Character string specifying the type of confidence interval.
#'   One of "wald", "percentile", "simul", or "none"
#' @param alpha Numeric significance level for confidence intervals (e.g., 0.05 for 95% CIs).
#' @param colors Character vector of colors. Length should match the number of
#'   unique methods in `plot_data`. If `NULL` or length doesn't match, ggplot2's
#'   default colors are used.
#'
#' @return A ggplot2 object with three faceted panels.
#'
#' @details
#' This is an internal helper function used by `plot.vefit()` and `compare_ve_fits()`.
#' The function automatically detects whether it's plotting a single method or
#' comparing multiple methods based on the `method` column in `plot_data`.
#'
#' For cumulative incidence panels, y-axis limits are shared to facilitate comparison.
#' The VE panel uses free y-axis scaling.
#'
#' @keywords internal
#' @noRd

plot_ve_panel <- function(plot_data,
                                   confint_type,
                                   alpha,
                                   colors = NULL) {

  n_methods <- length(unique(plot_data$method))

  # Check if plotting CIs
  has_ci <- (confint_type != "none")

  if (has_ci) {
    lower <- paste0(confint_type, "_lower")
    upper <- paste0(confint_type, "_upper")

    # These should exist if we got here (validation happened upstream)
    stopifnot(c(lower, upper) %in% names(plot_data))

    ci_level <- paste0((1 - alpha) * 100, "%")
    ci_label <- switch(confint_type,
                       "wald" = "pointwise Wald",
                       "percentile" = "pointwise percentile",
                       "simul" = "simultaneous")
    caption_text <- paste("\nShaded areas represent", ci_level, ci_label,
                          "confidence intervals.")
  } else {
    lower <- "estimate"
    upper <- "estimate"
    caption_text <- "\nPoint estimates only (no confidence intervals computed/requested)."
  }


  plot_data$term_label <- factor(plot_data$term, c("cuminc_0", "cuminc_1", "ve"),
                                 c("Cumulative Incidence:\n no vaccine",
                                   "Cumulative Incidence:\n vaccine",
                                   "VE"))

  # Calculate shared y-limits for cumulative incidence
  cuminc_data <- plot_data[plot_data$term %in% c("cuminc_0", "cuminc_1"), ]
  y_min <- min(cuminc_data[[lower]], na.rm = TRUE)
  y_max <- max(cuminc_data[[upper]], na.rm = TRUE)


  #Main plot
  p <- ggplot2::ggplot(plot_data,
                    ggplot2::aes(x = t0, y = estimate,
                                 color = method, fill = method)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~term_label, scales = "free") +
    ggh4x::facetted_pos_scales(
      y = list(
        ggplot2::scale_y_continuous(labels = scales::label_percent(),
                                    limits = c(y_min, y_max)),
        ggplot2::scale_y_continuous(labels = scales::label_percent(),
                                    limits = c(y_min, y_max)),
        ggplot2::scale_y_continuous(labels = scales::label_percent())
      )
    )

  # Add ribbon only if CIs exist
  if(has_ci){
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[lower]], ymax = .data[[upper]]),
                           alpha = if(n_methods > 1) 0.2 else 0.3,
                           linewidth = if(n_methods > 1) 0.1 else 0, show.legend = FALSE)
  }


  # Add color scales if appropriate, otherwise use default scales
  if (n_methods == length(colors)) {
    p <- p +
      ggplot2::scale_color_manual(values = colors, name = "Method") +
      ggplot2::scale_fill_manual(values = colors, name = "Method")
  }

  # Add labels and theme
  p <- p +
    ggplot2::labs(
      x = "Time since vaccination",
      y = "Estimate",
      caption = paste("\nShaded areas represent", ci_level, ci_label,
                      "confidence intervals.")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(hjust = 0),
      strip.background = ggplot2::element_rect(fill = "white", color = "white"),
      strip.text = ggplot2::element_text(size = 11, colour = "black"),
      axis.text = ggplot2::element_text(size = 9, colour = "black"),
      axis.title = ggplot2::element_text(size = 11, colour = "black")
    )

  # Add legend positioning if using method colors
  if (n_methods == 1) {
    p <- p + ggplot2::theme(
      legend.position = "none"
    )
  }else{
    p <- p + ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 10)
    )
  }

  return(p)
}


#' Validate and normalize CI type for plotting
#'
#' @param object A vefit object
#' @param confint_type Character string or NULL
#'
#' @return A validated confint_type string for plotting (one of "wald", "percentile", "simul")
#' @keywords internal
#' @noRd
validate_confint_type <- function(object, confint_type = NULL) {

  # Use object's confint_type if not specified
  if (is.null(confint_type)) {
    confint_type <- object$confint_type
  }

  # Handle "both" - default to "wald"
  if (confint_type == "both") {
    confint_type <- "wald"
  }

  # Validate it's one of the allowed types
  if (!confint_type %in% c("wald", "percentile", "simul", "none")) {
    stop("confint_type must be one of 'wald', 'percentile', 'simul', or 'none'")
  }

  # Check if simultaneous CIs are requested but not computed
  if (confint_type == "simul") {
    if (is.null(object$simul_z_star)) {
      stop("Simultaneous CIs not yet computed. Call add_simultaneous_ci()")
    }
  }

  if(confint_type == "none"){
    return("none")
  }else{
    #Check if requested CI type exists in estimates
    lower <- paste0(confint_type, "_lower")
    upper <- paste0(confint_type, "_upper")
    if (!all(c(lower, upper) %in% colnames(object$estimates$ve))) {
      stop("CI type '", confint_type, "' not found in object.\n")
    }
  }

  return(confint_type)
}



