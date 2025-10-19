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
#' Internal validation function that checks whether a user-supplied
#' list of marginalizing weights (`g_weights` and `p_weights`) has the required
#' structure and valid probability values.
#'
#' @param custom_weights List containing `g_weights` and `p_weights` data frames.
#' @param time_name Character; name of the vaccination time variable.
#' @param adjust_vars Character vector of covariate names.
#'
#' @return Invisibly returns `NULL` if all checks pass. Throws an error with a
#' descriptive message otherwise.
#' @keywords internal
#' @noRd
validate_marginalizing_weights <- function(custom_weights, time_name, adjust_vars) {
    # Check structure
    if (!is.list(custom_weights))
        stop("`custom_weights` must be a list containing 'g_weights' and 'p_weights'.", call. = FALSE)
    if (!all(c("g_weights", "p_weights") %in% names(custom_weights)))
        stop("List must include components 'g_weights' and 'p_weights'.", call. = FALSE)

    g <- custom_weights$g_weights
    p <- custom_weights$p_weights

    if (!is.data.frame(g) || !is.data.frame(p))
        stop("'g_weights' and 'p_weights' must be data frames.", call. = FALSE)

    # Required columns
    required_g <- c(time_name, "prob", adjust_vars)
    required_p <- c("prob", adjust_vars)

    check_missing_cols <- function(df, required, name) {
        missing <- setdiff(required, names(df))
        if (length(missing))
            stop(sprintf("`%s` missing required columns: %s",
                         name, paste(missing, collapse = ", ")), call. = FALSE)
    }

    check_missing_cols(g, required_g, "g_weights")
    check_missing_cols(p, required_p, "p_weights")

    # Validate probabilities
    validate_prob_column(g$prob, "g_weights")
    validate_prob_column(p$prob, "p_weights")

    # Check sums
    if (any(abs(tapply(g$prob, g[adjust_vars], sum) - 1) > 1e-6))
        stop("In 'g_weights', probabilities must sum to 1 within each group.", call. = FALSE)

    if (abs(sum(p$prob) - 1) > 1e-6)
        stop("In 'p_weights', probabilities must sum to 1.", call. = FALSE)

    invisible(NULL)
}

#' Validate numeric probability vector
#' @keywords internal
#' @noRd
validate_prob_column <- function(prob, name) {
    if (!is.numeric(prob))
        stop(sprintf("'prob' column in '%s' must be numeric.", name), call. = FALSE)
    if (any(is.na(prob)))
        stop(sprintf("'prob' column in '%s' contains missing values.", name), call. = FALSE)
    if (any(prob < 0 | prob > 1))
        stop(sprintf("'prob' column in '%s' must be between 0 and 1.", name), call. = FALSE)
    invisible(NULL)
}

# g <- weights$g_weights
# g$prob[g$x2 == 5] <- 0
# p_weights <- weights$p_weights
# p_weights$prob[p_weights$x2 == 5] <- NA
# p <- p_weights

canonicalize_weights <- function(weights, time_name, adjust_vars) {
    g <- weights$g_weights
    p <- weights$p_weights

    # Rename prob columns if needed
    if (!"prob_g" %in% names(g)) names(g)[names(g)=="prob"] <- "prob_g"
    if (!"prob_p" %in% names(p)) names(p)[names(p)=="prob"] <- "prob_p"

    # Drop rows with non-positive weight
    gp <- merge(g, p, all = TRUE)
    prob_g_zero <- gp$prob_g == 0 & !is.na(gp$prob_g)
    prob_p_zero <- gp$prob_p == 0 & !is.na(gp$prob_p)

    gp <- gp[!prob_g_zero & !prob_p_zero,]

    if(anyNA(gp$prob_g * gp$prob_p)){
        stop("Covariate groups differ between g_weights and p_weights.", call. = FALSE)
    }
    # Add group_id
    groups <- unique(gp[adjust_vars])
    groups <- groups[do.call(order, groups[adjust_vars]), , drop = FALSE]
    groups$group_id <- seq_len(nrow(groups))

    g <- merge(groups, g, all.x = TRUE)
    p <- merge(groups, p, all.x = TRUE)

    # Return canonical structure
    list(
        g_weights = g[, c("group_id", adjust_vars, time_name, "prob_g"), drop = FALSE],
        p_weights = p[, c("group_id", adjust_vars, "prob_p"), drop = FALSE]
    )
}



#' @keywords internal
check_reserved_vars <- function(vars, reserved_vars, vars_label) {
    conflicts <- intersect(reserved_vars, vars)
    if (length(conflicts) > 0) {
        stop(
            paste(vars_label, " contain reserved variable names: "),
            paste(conflicts, collapse = ", "),
            "\nPlease rename these variables before fitting the model.",
            call. = FALSE
        )
    }
}

