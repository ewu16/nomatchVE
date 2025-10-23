#' Match vaccinated and unvaccinated individuals using rolling cohort design
#'
#' @description Creates 1:1 matched pairs of vaccinated ("cases") and
#' unvaccinated ("controls") individuals. Uses a rolling cohort design where
#' controls must be unvaccinated and event-free at the time they are matched to
#' a case.
#'
#' @details For each vaccination time, newly vaccinated individuals are matched
#' to eligible controls using exact covariate matching. Controls are eligible if
#' unvaccinated and event-free at that time. Vaccinated individuals may appear
#' as a control (when they are not yet vaccinated) and as a case.
#'
#'
#' @inheritParams nomatchVE
#' @param data Data frame with study population
#' @param id_name Name of unique identifier variable of individuals
#' @param matching_vars Character vector of variables to match on exactly
#' @param replace Logical. Allow controls to be reused? Default: `FALSE`
#' @param seed Integer for reproducibility. Default: `NULL`


#' @return A list containing the following:
#' \describe{
#'   \item{matched_data}{Data frame of matched pairs with original variables plus:
#'       \itemize{
#'       \item `match_index_time`: Matching time
#'       \item `match_type`: "case" or "control"
#'       \item `match_<exposure>`: Treatment at matching
#'       \item `match_id`: Pair identifier
#'     }
#'    }
#'   \item{n_unmatched_cases}{Number of unmatched vaccinated individuals}
#'   \item{discarded}{Logical vector indicating excluded individuals}
#' }
#'
#' @export
#' @examples
#' matched_cohort <- match_rolling_cohort(
#' data = simdata,
#' outcome_time =  "Y",
#' exposure = "V",
#' exposure_time = "D_obs",
#' matching_vars = c("x1", "x2"),
#' id_name = "ID",
#' seed = 5678
#' )
match_rolling_cohort <- function(data, outcome_time, exposure, exposure_time, matching_vars, id_name, replace = FALSE, seed = NULL){

    # Check data/inputs
    validate_data(
        data = data,
        outcome_time = outcome_time,
        exposure = exposure,
        exposure_time = exposure_time,
        covariates = matching_vars
        )

    if(any(duplicated(data[[id_name]]))){
        stop("<id_name> is not a unique identifier for data")
    }

    # Check for name conflicts for variables that will be added
    reserved_vars <- c("match_index_time", "match_type", paste0("match_", exposure), "match_id")

    check_reserved_vars(names(data), reserved_vars, "Names in data")


    #set seed for reproducibility due to random matching
    if(!is.null(seed)){
        set.seed(seed)
    }
    #store information about matching procedure
    used_controls <- vector() #store IDs of controls used
    matched_df_list <- list() #store matched groups
    n_pairs_counter <- 0 #keep track of how many pairs have been matched

    #vaccinated serve as "cases" for matching, unvaccinated are "controls"
    cases <- subset(data, data[[exposure]] == 1)
    observed_ds <- sort(unique(data[[exposure_time]]))

    #loop through observed vaccination eval_times- each newly vaccinated person is eligible for matching
    for(d in observed_ds){
        #cases on day d
        case_d <- subset(cases, cases[[exposure_time]] == d)

        #controls on day d
        ## eligibility conditions
        unvaxed_on_d <- data[[exposure_time]] > d | is.na(data[[exposure_time]])
        endpointfree_on_d <- data[[outcome_time]] > d
        eligible_on_d <- unvaxed_on_d & endpointfree_on_d
        control_d <- subset(data, eligible_on_d)

        ##further restrict control pool if resampling of controls not allowed
        if(replace == FALSE){
           control_d <- subset(control_d, !control_d[[id_name]] %in% used_controls)
        }

        #get matches for day d within each subgroup defined by matching variables
        case_d_split <- split(case_d, case_d[,matching_vars], drop = TRUE)
        control_d_split <- split(control_d, control_d[,matching_vars], drop = TRUE)
        #subgroups that exist for both cases and controls
        subgroup_names <- intersect(names(case_d_split), names(control_d_split))

        #loop through subgroups
        for(subgroup_name in subgroup_names){
            case_d_subgroup <- case_d_split[[subgroup_name]]
            control_d_subgroup <- control_d_split[[subgroup_name]]

            #choose cases and controls
            n_matches <- min(nrow(case_d_subgroup),  nrow(control_d_subgroup))
            case_inds <- 1:n_matches
            #case_inds <- sample(1:nrow(case_d_subgroup), size = n_matches, replace = FALSE)
            control_inds <- sample(1:nrow(control_d_subgroup), size = n_matches, replace = FALSE)

            #stack cases and controls to create matched df
            matched_cases <- case_d_subgroup[case_inds,]
            matched_controls <- control_d_subgroup[control_inds,]
            matched_df <- rbind(matched_cases, matched_controls)

            #add information about matched pairs to matched df
            matched_df$match_index_time <- d
            matched_df$match_type <- rep(c("case", "control"), each = n_matches)
            matched_df[[paste0("match_", exposure)]] <- rep(c(1, 0), each = n_matches)
            matched_df$match_id <- rep(n_pairs_counter + (1:n_matches), 2)

            #update info for larger for loops
            matched_df_list <- c(matched_df_list, list(matched_df))
            used_controls <- c(used_controls, matched_controls[[id_name]])
            n_pairs_counter <- n_pairs_counter + n_matches

        }#end loop for covariate groups
    }#end loop for observed vaccination eval_times

    #handle case where no matches are found
    if(length(matched_df_list) == 0) {
        warning("No matches found. Consider relaxing matching criteria.")
        return(list(matched_data = data.frame(),
                    n_unmatched_cases = nrow(cases),
                    discarded = rep(TRUE, nrow(data))))
    }

    #otherwise, compile final matched dataset
    matched <- do.call(rbind, matched_df_list)
    matched_final <- matched[order(matched$match_id),]
    rownames(matched_final) <- NULL

    unmatched_cases <- setdiff(cases[[id_name]], matched_final[matched_final$match_type == "case",][[id_name]])
    discarded <- !(data[[id_name]] %in% matched_final[[id_name]])

    list(matched_data = matched_final,
         n_unmatched_cases = length(unmatched_cases),
         discarded = discarded)
}
