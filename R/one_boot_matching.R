#'  Compute one bootstrap replicate of matching-based VE point estimate
#'
#' @inheritParams get_one_matching_ve
#'
#' @return  A matrix of bootstrapped estimates where the the columns of the matrix are the cumulative
#'  incidence/VE terms and the rows are the requested time points for evaluation.
#'
#' @keywords internal
#'
one_boot_matching <- function(matched_data,
                              outcome_time,
                              outcome_status,
                              exposure,
                              exposure_time,
                              tau,
                              eval_times,
                              limit_type = "fixed",
                              data = NULL,
                              id_name = NULL,
                              matching_vars = NULL,
                              replace = NULL,
                              keep_boot_samples = TRUE){



    # --------------------------------------------------------------------------
    # 1. Create bootstrapped sample(s)
    # --------------------------------------------------------------------------
    if(limit_type == "fixed"){
        #Bootstrap from fixed matched cohort
        boot_matched_ids <- sample(unique(matched_data$match_id), replace = TRUE)
        boot_matched_inds <- as.vector(sapply(boot_matched_ids, \(match_id) which(matched_data$match_id == match_id)))
        boot_matched_data <- matched_data[boot_matched_inds,]
        boot_matched_data$match_id <- rep(1:length(boot_matched_ids), each = 2)

    }else if(limit_type == "limit"){
        #Check required args are provided
        required_args <- list(data = data,
                              id_name = id_name,
                              matching_vars = matching_vars,
                              replace = replace)

        missing_args <- names(required_args)[vapply(required_args, is.null, logical(1))]
        if (length(missing_args) > 0) {
            stop(sprintf("When limit_type = 'limit', the following arguments must be specified: %s",
                         paste(missing_args, collapse = ", "))
            )
        }

        #Bootstrap the original data and rematch
        boot_inds <- sample(1:nrow(data), replace = TRUE)
        boot_data <-  data[boot_inds,]
        boot_id_name <- paste0("boot_",id_name)
        boot_data[[boot_id_name]] <- 1:nrow(boot_data)

        boot_matched_cohort <- match_rolling_cohort(
            data = boot_data,
            outcome_time = outcome_time,
            exposure = exposure,
            exposure_time = exposure_time,
            id_name =  boot_id_name,
            matching_vars = matching_vars,
            replace = replace
            )

        boot_matched_data <- boot_matched_cohort[[1]]

    }

    # --------------------------------------------------------------------------
    # 2. Compute VE for bootstrapped data
    # --------------------------------------------------------------------------
    boot_matching_ve <- get_one_matching_ve(matched_data = boot_matched_data,
                                            outcome_time = outcome_time,
                                            outcome_status = outcome_status,
                                            exposure = exposure,
                                            exposure_time = exposure_time,
                                            tau = tau,
                                            eval_times = eval_times,
                                            keep_models = FALSE,
                                            keep_adata = FALSE)


    boot_matching_ve$pt_estimates
}




