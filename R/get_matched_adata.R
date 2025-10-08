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
#' @param outcome_name Character string representing the original outcome variable in  `matched_data`
#' @param event_name Character string representing the original event variable in `matched_data`
#' @param trt_name Character string representing the original vaccination status variable in `matched_data`
#' @param time_name Character string representing the original vaccination time variable in `matched_data`
#' @param tau The time excluded after vaccination to allow building up of
#'   immunity
#' @param pair_censoring Logical to indicate if matched pairs should be censored
#'   when unvaccinated pair becomes vaccinated.
#' @return A data frame representing the analysis-eligible matched dataset (a subset of `matched data`).
#' Contains all variables in `matched_data`, plus the following variables
#' \describe{
#' \item{<event_name>_d}{Event variable after adjusting for pair_censoring}
#' \item{<outcome_name>_matched}{Outcome variable after adjusting for pair_censoring}
#' \item{T_d}{Time to endpoint from matched index date (adjusted for pair_censoring if requested)}
#' }
#' @export
#'
clean_matched_data <- function(matched_data, outcome_name, event_name, trt_name, time_name,
                               tau, pair_censoring = TRUE){
    if(pair_censoring == FALSE){
        tmp_data <- matched_data
        y_name <- outcome_name
    }else{
        #Censor cases and controls if control later gets vaccinated
        ordered_data <- matched_data[order(matched_data$pair_id),]
        #new variables
        new_outcome_name <- paste0(outcome_name, "_matched")
        new_event_name <- paste0(event_name, "_d")
        ordered_data[[new_outcome_name]] <-ordered_data[[new_event_name]] <-  NA

        cases <- subset(ordered_data, type == "case")
        controls <- subset(ordered_data, type == "control")
        stopifnot(all.equal(cases$pair_id, controls$pair_id))

        #controls who later got vaccinated
        trt_later <- controls[[trt_name]] == 1

        #censor controls at their time of vaccination
        controls[[new_outcome_name]] <- ifelse(trt_later, controls[[time_name]], controls[[outcome_name]])
        controls[[new_event_name]]   <- ifelse(trt_later, 0, controls[[event_name]])

        #censor cases at control's time of vaccination if still at risk
        at_risk <-   cases[[outcome_name]] >  controls[[time_name]]
        cases[[new_outcome_name]] <- ifelse(trt_later & at_risk, controls[[time_name]], cases[[outcome_name]])
        cases[[new_event_name]]   <- ifelse(trt_later & at_risk, 0, cases[[event_name]])

        tmp_data <- rbind(cases, controls)
        y_name <- new_outcome_name
    }

    #Create matching analysis dataset
    ##compute survival time from matching index date
    tmp_data$T_d <- tmp_data[[y_name]] - tmp_data$d
    ##exclude pairs where at least one individual is no longer at risk by tau
    excluded_pair_ids <- unique(tmp_data$pair_id[!tmp_data$T_d > tau])
    adata <- subset(tmp_data, !tmp_data$pair_id %in% excluded_pair_ids)
    adata[order(adata$pair_id),]

}


