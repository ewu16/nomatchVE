# Format estimates --------------------------------------------------------

#' Convert estimates to tidy long-format data frame
#'
#' @description
#' Converts the `estimates` component of a `vefit` object into a
#' long-format data frame for easy plotting and analysis.
#' Each row represents one term (cumulative incidence or derived effect measure) at a specific timepoint.
#'
#' @param estimates A list of matrices, typically from `fit$estimates` in a `vefit`
#'   object. Each matrix (e.g., `cuminc_0`, `cuminc_1`, `vaccine_effectiveness`/`risk_ratio`)
#'   contains estimates, confidence intervals, and related statistics evaluated
#'   at different timepoints.
#'
#' @return
#' A long-format data frame with:
#'   - `t0`: the timepoint
#'   - `term`: what was estimated (`cuminc_0`, `cuminc_1`, or `vaccine_effectiveness`/`risk_ratio`)
#'   -  Remaining columns from the original matrices (estimates, confidence intervals, ad related statistics.)
#'
#' @export
#'
estimates_to_df <- function(estimates){
    df_list <- lapply(seq_along(estimates), \(i){
        df <-data.frame(estimates[[i]])
        df$term <- names(estimates)[i]
        df$t0 <- as.numeric(rownames(df))
        df <- df[, c("t0", "term", colnames(estimates[[i]]))]
    })
    estimates_df <- do.call("rbind", df_list)
    ordered_df <-  estimates_df[order(estimates_df$t0),]
    rownames(ordered_df) <- NULL
    ordered_df

}

