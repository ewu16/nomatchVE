#' Convert list of estimates to dataframe with t0 as column
#'
#' @param object Object created by timepoints
#'
#' @return Data frame with estimates across all time points
#' @export
#'
timepoints_to_df <- function(object){
  estimates_list <- lapply(seq_along(object$estimates), \(i){
    estimates_df <- estimates_to_df(object$estimates[[i]])
    estimates_df$t0 <- names(object$estimates)[i]
    estimates_df
  })
  do.call("rbind", estimates_list)
}


#' Plot cumulative incidence/VE estimates
#'
#' @param object Object created by timepoints
#'
#' @return ggplot
#' @export
plot_timepoints <- function(object){
  #TODO:specify what confidence interval types to plot
  df <- timepoints_to_df(object)
  df |>
    ggplot2::ggplot(aes(x = t0, y = estimate,group = 1)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = wald_lower, ymax = wald_upper),
                alpha = 0.15, linewidth = 0, show.legend = FALSE) +
    ggplot2::facet_wrap(~term, scales = "free") +
    ggplot2::scale_y_continuous(labels = scales::label_percent()) +
    ggplot2::labs(x = "Time since vaccination (t)", y = "",
         color = "Method:") +
    ggplot2:: theme_bw() +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0.5),
                   strip.background = ggplot2::element_rect(fill = "white",color = "white"),
                   strip.text = ggplot2::element_text(size = 12))

}
