# Format estimates --------------------------------------------------------


#' Convert estimates list to tidy data frame
#'
#' @param estimates {A list of matrices for cumulative incidences and VE.
#'  Each matrix contains the point estimate and confidence intervals for the specified term.}
#'
#' @return A data frame with columns `t0` specifying the time point and `term`,
#' and point estimates, and confidence intervals.
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


#' Reform estimates list to a list arranged by time point rather than term
#'
#' @param estimates {A list of matrices for cumulative incidences and VE.
#'  Each matrix contains the point estimate and confidence intervals for the specified term.}
#'
#'
#' @return A list of matrices for each time point where the rows
#' are the term and the columns are the estimates at the given time point.
#' @export
estimates_by_time <- function(estimates){
    times <- rownames(estimates[[1]])
    estimates_time_list <- lapply(seq_along(times), \(i){
        tmp <- rbind(estimates[[1]][i,],
                     estimates[[2]][i,],
                     estimates[[3]][i,])
        rownames(tmp) <- names(estimates)
        tmp
    })
    names(estimates_time_list) <- times
    estimates_time_list
}




