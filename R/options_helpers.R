#' Parameters for creating a matched dataset and matched analysis data set
#'
#' @inheritParams match_rolling_cohort
#' @inheritParams clean_matched_data
#'
#' @return A list with control parameters
#' @keywords internal
#'
matched_dist <- function(id_name = "ID",
                         matched_data = NULL,
                         replace = FALSE,
                         seed = NULL){
    list(id_name = id_name,
         matched_data = matched_data,
         replace = replace,
         seed = seed)
}

