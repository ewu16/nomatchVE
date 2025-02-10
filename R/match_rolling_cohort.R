#' Perform exact 1:1 matching using "rolling cohort" design
#'
#' @description This function creates a matched cohort based on the rolling cohort design.
#'    Individuals are eligible to be selected as unvaccinated "controls" only if they
#'    are unvaccinated and endpoint free on the potential matching date
#'
#' @inheritParams obsve
#' @param data A data frame for which matches are sought
#' @param id_name Character string specifying name of a unique individual
#'   identifier in `data`.
#' @param matching_vars A character vector containing the names of variables in
#'   `data` to match on .
#' @param replace Logical: Should matching be done with replacement?
#' @param seed Numeric: seed to set prior to performing random matching



#' @return A list containing the following:
#' \describe{
#' \item{matched_data}{A data frame containing the matched individuals}
#' \item{n_unmatched_cases}{The number of vaccinated individuals for whom no match was found}
#' \item{discarded}{A logical vector denoting whether individual in the source population was
#' not included in matched_cohort}
#' }
#' @export
#'
match_rolling_cohort <- function(data, outcome_name, trt_name, time_name, id_name, matching_vars, replace = FALSE, seed = NULL){
    #set seed for reproducibility due to random matching
    if(!is.null(seed)){
        set.seed(seed)
    }
    #store information about matching procedure
    used_controls <- vector() #store IDs of controls used
    matched_df_list <- list() #store matched groups
    n_pairs_counter <- 0 #keep track of how many pairs have been matched

    #vaccinated serve as "cases" for matching, unvaccinated are "controls"
    cases <- subset(data, data[[trt_name]] == 1)
    observed_ds <- sort(unique(data[[time_name]]))

    #loop through observed vaccination times- each newly vaccinated person is eligible for matching
    for(d in observed_ds){
        #cases on day d
        case_d <- subset(cases, cases[[time_name]] == d)

        #controls on day d
        ## eligibility conditions
        unvaxed_on_d <- data[[time_name]] > d | is.na(data[[time_name]])
        endpointfree_on_d <- data[[outcome_name]] > d
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
            matched_df$d <- d
            matched_df$type <- rep(c("case", "control"), each = n_matches)
            matched_df[[paste0(trt_name, "_d")]] <- rep(c(1, 0), each = n_matches)
            matched_df$pair_id <- rep(n_pairs_counter + (1:n_matches), 2)

            #update info for larger for loops
            matched_df_list <- c(matched_df_list, list(matched_df))
            used_controls <- c(used_controls, matched_controls[[id_name]])
            n_pairs_counter <- n_pairs_counter + n_matches

        }#end loop for covariate groups
    }#end loop for observed vaccination times

    #compile final matched dataset
    matched <- do.call(rbind, matched_df_list)
    matched_final <- matched[order(matched$pair_id),]
    rownames(matched_final) <- NULL

    unmatched_cases <- setdiff(cases[[id_name]], subset(matched_final, type == "case")[[id_name]])
    discarded <- !(data[[id_name]] %in% matched_final[[id_name]])

    list(matched_data = matched_final,
         n_unmatched_cases = length(unmatched_cases),
         discarded = discarded)
}
