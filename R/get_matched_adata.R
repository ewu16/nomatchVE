#' Get analysis matched dataset from matched cohort
#'
#' @description This function modifies the original matched data, preparing it for use in
#' analysis. Namely, this includes
#' *  (if requested) censoring matched pairs at the time
#' of the unvaccinated individual's vaccination
#' * creating the outcome time from matching index date to endpoint
#' * restricting to matched  pairs in which both individuals are eligible `tau` days after the matching
#' index date
#'

#' @param matched_data A data frame representing a matched cohort
#' @param outcome_time Character string representing the original outcome variable in  `matched_data`
#' @param outcome_status Character string representing the original event variable in `matched_data`
#' @param exposure Character string representing the original vaccination status variable in `matched_data`
#' @param exposure_time Character string representing the original vaccination time variable in `matched_data`
#' @param tau The time excluded after vaccination to allow building up of
#'   immunity
#' @return A data frame representing the analysis-eligible matched dataset (a subset of `matched data`).
#' Contains all variables in `matched_data`, plus the following variables
#' \describe{
#' \item{match_<outcome_status>}{Event variable after adjusting for pair_censoring}
#' \item{match_<outcome_time>}{Outcome variable after adjusting for pair_censoring}
#' \item{match_T}{Time to endpoint from matched index date after adjusting for pair_censoring}
#' }
#' @keywords internal
#'
clean_matched_data <- function(matched_data, outcome_time, outcome_status, exposure, exposure_time,
                               tau){

    # Check for name conflicts for variables that will be added
    reserved_vars <- c("match_T", paste0("match_", outcome_time), paste0("match_", outcome_status))
    conflicts <- intersect(reserved_vars, names(matched_data))

    if(length(conflicts) > 0) {
        stop(
            "Data contains reserved variable names: ",
            paste(conflicts, collapse = ", "), "\n",
            "Please rename these variables before matching.",
            call. = FALSE
        )
    }

    #Censor cases and controls if control later gets vaccinated
    ordered_data <- matched_data[order(matched_data$match_id),]
    #new variables
    new_outcome_name <- paste0("match_", outcome_time)
    new_event_name <- paste0("match_", outcome_status)
    ordered_data[[new_outcome_name]] <-ordered_data[[new_event_name]] <-  NA

    cases <- ordered_data[ordered_data$match_type == "case",]
    controls <-  ordered_data[ordered_data$match_type == "control",]
    stopifnot(all.equal(cases$match_id, controls$match_id))

    #controls who later got vaccinated
    trt_later <- controls[[exposure]] == 1

    #censor controls at their time of vaccination
    controls[[new_outcome_name]] <- ifelse(trt_later, controls[[exposure_time]], controls[[outcome_time]])
    controls[[new_event_name]]   <- ifelse(trt_later, 0, controls[[outcome_status]])

    #censor cases at control's time of vaccination if still at risk
    at_risk <-   cases[[outcome_time]] >  controls[[exposure_time]]
    cases[[new_outcome_name]] <- ifelse(trt_later & at_risk, controls[[exposure_time]], cases[[outcome_time]])
    cases[[new_event_name]]   <- ifelse(trt_later & at_risk, 0, cases[[outcome_status]])

    tmp_data <- rbind(cases, controls)

    #Create matching analysis dataset
    ##compute survival time from matching index date
    tmp_data$match_T <- tmp_data[[new_outcome_name]] - tmp_data$match_index_time
    ##exclude pairs where at least one individual is no longer at risk by tau
    excluded_match_ids <- unique(tmp_data$match_id[!tmp_data$match_T > tau])
    adata <- subset(tmp_data, !tmp_data$match_id %in% excluded_match_ids)

    adata[order(adata$match_id),]

}


