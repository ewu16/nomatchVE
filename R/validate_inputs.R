#' Validate inputs to main functions
#'
#' @description
#' Internal validation function to check that all inputs to `nomatchVE()/matching_ve()` are
#' properly formatted and logically consistent.
#'
#' @inheritParams nomatchVE
#'
#' @return Invisibly returns `NULL` if all checks pass. Throws an error with
#'   descriptive message if any validation fails.
#'
#' @keywords internal
#' @noRd
validate_ve_inputs <- function(data, outcome_name, event_name, trt_name, time_name, adjust_vars, times, tau, censor_time) {

    validate_data(data, outcome_name, event_name, trt_name, time_name, adjust_vars)


    # Check time relationships
    if(tau < 0){
        stop("Tau must be non-negative (>= 0)")
    }
    if(tau > 28) {
        warning(
            "tau = ", tau, " days is unusually long for an immune build-up period.\n",
            "Resulting estimates may require strong assumptions for causal interpretation.\n",
            "Please verify this is intentional."
        )
    }

    if(any(times <= tau)) {
        stop("All evaluation times 'times' must be greater than tau (", tau, " days)")
    }
    if(length(censor_time) != 1){
        stop("Censoring time 'censor_time' must be a scalar value")
    }
    if(censor_time < max(times)){
        stop("Censoring time 'censor_time' must be greater than or equal to max(times)")
    }
}

validate_data <- function(data, outcome_name, event_name, trt_name, time_name, adjust_vars){
    # Check data is a data frame
    if(!is.data.frame(data)) {
        stop("'data' must be a data.frame, not ", class(data)[1])
    }
    if(nrow(data) == 0) {
        stop("'data' has no rows")
    }

    # Check required columns exist
    required_cols <- c(outcome_name, event_name, trt_name, time_name, adjust_vars)
    missing_cols <- setdiff(required_cols, names(data))
    if(length(missing_cols) > 0) {
        stop("Missing required column(s) in data: ", paste(missing_cols, collapse = ", "))
    }

    # Check missing values
    if(sum(is.na(data[, c(outcome_name, event_name, trt_name)])) > 0){
        stop("Missing values in outcome, event or treatment variables are not supported. Please remove these observations and try again.")
    }
    if(sum(is.na(data[, adjust_vars])) > 0){
        stop("Missing values in adjustment variables are not supported. Please remove these observations and try again.")
    }

    # Check data types and coding
    if(!all(data[[event_name]] %in% c(0, 1))) {
        stop("Event variable '", event_name, "' must be coded as 0/1 (0=censored, 1=event)")
    }
    if(!all(data[[trt_name]] %in% c(0, 1))) {
        stop("Treatment variable '", trt_name, "' must be coded as 0/1 (0=unvaccinated, 1=vaccinated)")
    }
    if(any(data[[outcome_name]] < 0)){
        stop("Outcome variable '", outcome_name, "' must be non-negative")
    }
    if(any(data[[time_name]][!is.na(data[[time_name]])] < 0)){
        stop("Time of treatment variable '", time_name, "' must be non-negative")
    }

    #Check vaccination times
    vaccinated <- data[[trt_name]] == 1
    vax_time <- data[[time_name]]
    if(any(is.na(vax_time[vaccinated]))){
        stop("All vaccinated individuals must have a non-missing vaccination time")
    }
    if(any(!is.na(vax_time[!vaccinated]))){
        stop("All unvaccinated individuals must have a missing (NA) vaccination time")
    }
}


#' Validate custom marginalizing weights
#'
#' @description
#' Internal validation function to check that user-supplied marginalizing
#' weights (`g_dist` and `p_dist`) are properly formatted with required
#' columns and valid probability values.
#'
#' @param custom_weights Either a character string or a list containing
#'   `g_dist` and `p_dist` data frames. Only list inputs are validated by
#'   this function.
#' @param time_name Character string specifying the vaccination time variable name
#' @param adjust_vars Character vector of covariate names
#'
#' @return Invisibly returns `NULL` if all checks pass. Throws an error with
#'   descriptive message if any validation fails.
#'
#' @keywords internal
#' @noRd
validate_marginalizing_weights <- function(custom_weights, time_name, adjust_vars) {

    # Only validate if it's a list (custom weights)
    # Character values are validated by match.arg() in main function
    if(!is.list(custom_weights)) {
        stop("Custom marginalizing weights must be a list")
    }

    # Check list structure
    if(!all(c("g_dist", "p_dist") %in% names(custom_weights))) {
        stop(
            "Custom marginalizing weights must be a list with 'g_dist' and 'p_dist'.\n",
            "Found components: ", paste(names(custom_weights), collapse = ", "),
            call. = FALSE
        )
    }

    g_dist <- custom_weights$g_dist
    p_dist <- custom_weights$p_dist

    # Validate structure
    if(!is.data.frame(g_dist) || !is.data.frame(p_dist)) {
        stop("Both g_dist and p_dist must be data frames", call. = FALSE)
    }

    # Check required columns
    required_g <- c("group_name", time_name, "prob", adjust_vars)
    missing_g <- setdiff(required_g, names(g_dist))
    if(length(missing_g) > 0) {
        stop(
            "g_dist missing columns: ", paste(missing_g, collapse = ", "),
            call. = FALSE
        )
    }

    required_p <- c("group_name", "prob", adjust_vars)
    missing_p <- setdiff(required_p, names(p_dist))
    if(length(missing_p) > 0) {
        stop(
            "p_dist missing columns: ", paste(missing_p, collapse = ", "),
            call. = FALSE
        )
    }

    # Validate probabilities
    validate_g_dist_probabilities(g_dist)
    validate_p_dist_probabilities(p_dist)

    invisible(NULL)
}

#' Validate g_dist probabilities sum to 1 within each covariate group
#'
#' @param g_dist Data frame with columns `group_name`, `prob`, and covariate columns
#'
#' @return Invisibly returns `NULL` if all checks pass. Throws an error if
#'   probabilities don't sum to 1 within groups or contain invalid values.
#'
#' @keywords internal
#' @noRd
validate_g_dist_probabilities <- function(g_dist) {

    # Check probabilities are numeric
    if(!is.numeric(g_dist$prob)) {
        stop("'prob' column in 'g_dist' must be numeric", call. = FALSE)
    }

    # Check probabilities are between 0 and 1
    if(any(g_dist$prob < 0 | g_dist$prob > 1, na.rm = TRUE)) {
        stop("'prob' in 'g_dist' must be between 0 and 1", call. = FALSE)
    }

    # Check for missing probabilities
    if(any(is.na(g_dist$prob))) {
        stop("'prob' in 'g_dist' contains missing values", call. = FALSE)
    }

    # Calculate sum of probabilities within each group
    prob_sums <- tapply(g_dist$prob, g_dist$group_name, sum)

    # Check if sums are approximately 1 (allow small floating point errors)
    tolerance <- 1e-6
    bad_groups <- names(prob_sums)[abs(prob_sums - 1) > tolerance]

    if(length(bad_groups) > 0) {
        # Show details for first few problematic groups
        n_show <- min(3, length(bad_groups))
        examples <- paste0(
            bad_groups[1:n_show],
            " (sum=", round(prob_sums[bad_groups[1:n_show]], 4), ")",
            collapse = ", "
        )
        if(length(bad_groups) > n_show) {
            examples <- paste0(examples, ", ...")
        }

        stop(
            "In 'g_dist', probabilities must sum to 1 within each covariate group.\n",
            "Groups with incorrect sums: ", examples, "\n",
            "Total problematic groups: ", length(bad_groups),
            call. = FALSE
        )
    }

    invisible(NULL)
}


#' Validate p_dist probabilities sum to 1
#'
#' @param p_dist Data frame with columns `group_name`, `prob`, and covariate columns
#'
#' @return Invisibly returns `NULL` if all checks pass. Throws an error if
#'   probabilities don't sum to 1 or contain invalid values.
#'
#' @keywords internal
#' @noRd
validate_p_dist_probabilities <- function(p_dist) {

    # Check probabilities are numeric
    if(!is.numeric(p_dist$prob)) {
        stop("'prob' column in 'p_dist' must be numeric", call. = FALSE)
    }

    # Check probabilities are between 0 and 1
    if(any(p_dist$prob < 0 | p_dist$prob > 1, na.rm = TRUE)) {
        stop("'prob' in 'p_dist' must be between 0 and 1", call. = FALSE)
    }

    # Check for missing probabilities
    if(any(is.na(p_dist$prob))) {
        stop("'prob' in 'p_dist' contains missing values", call. = FALSE)
    }

    # Check total sum
    total_prob <- sum(p_dist$prob)
    tolerance <- 1e-6

    if(abs(total_prob - 1) > tolerance) {
        stop(
            "In 'p_dist', probabilities must sum to 1.\n",
            "Current sum: ", round(total_prob, 6),
            call. = FALSE
        )
    }

    invisible(NULL)
}
