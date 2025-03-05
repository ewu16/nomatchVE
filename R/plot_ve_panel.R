#' Plot cumulative incidence/VE estimates
#'
#' @param plot_data A dataframe
#' @param ci_type String indicating type of confidence interval to plot.
#' Takes values of  "wald", "percentile" or "simul".
#'
#' @return ggplot
#' @export
plot_ve_panel <- function(object, ci_type){
  lower <- paste0(ci_type, "_lower")
  upper <- paste0(ci_type, "_upper")


  plot_data <- estimates_to_df(object$estimates)
  stopifnot(c(lower, upper) %in% names(plot_data))

  plot_data$term_label <- factor(plot_data$term, c("risk_0", "risk_1", "ve"),
                                 c("Cumulative incidence: no vaccine",
                                   "Cumulative incidence: vaccine",
                                   "VE"))


  plot_data |>
    ggplot2::ggplot(ggplot2::aes(x = t0, y = estimate,group = 1)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[lower]], ymax = .data[[upper]]),
                alpha = 0.5, linewidth = 0, show.legend = FALSE) +
    ggplot2::facet_wrap(~term_label, scales = "free") +
    ggplot2::scale_y_continuous(labels = scales::label_percent()) +
    ggplot2::labs(x = "Time since vaccination", y = "Estimate") +
    ggplot2:: theme_bw() +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0.5),
                   strip.background = ggplot2::element_rect(fill = "white",color = "white"),
                   strip.text = ggplot2::element_text(size = 10))
}
