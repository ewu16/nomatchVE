#' Estimate observed distributions of exposure times and covariates
#'
#' @description Computes two empirical probability distributions used to marginalize
#' time- and covariate- specific cumulative incidences:
#' - \eqn{g(d \mid x)}: the distribution of exposure times within each covariate subgroup
#' of exposed individuals who remain at risk `tau` days after exposure.
#' - \eqn{p(x)}: the distribution of covariates among exposed individuals
#' who remain at risk `tau` days after exposure.
#'
#' Called internally by [nomatchVE()]. Provided as an example of the
#' structure needed if a user passes in `custom_weights`.
#'
#' @inheritParams get_gp
#'
#' @return A list with two data frames:
#' - `g_weights`: covariate-conditional exposure time  probabilities (\eqn{g(d \mid x)})
#' - `p_weights`: covariate probabilities (\eqn{p(x)})
#'
#'   Each data frame includes all variables in `adjust_vars`, and a `prob` column
#'   with empirical probabilities.
#'   `g_weights` additionally includes a `<time_name>` column for exposure times.
#'
#'
#' @export
#' @examples
#' weights <- get_observed_weights(simdata, "Y", "V", "D_obs",
#'                                    c("x1", "x2"), tau = 14)
#' str(weights)
#'
get_observed_weights <- function(data, outcome_name, trt_name,
                            time_name, adjust_vars, tau){

    gp_list <- get_gp(data, outcome_name, trt_name, time_name, adjust_vars, tau)

    #Simplify output to be similar to user-provided weights
    g <- gp_list$g_weights
    p <- gp_list$p_weights

    # Drop group_id
    g <- g[, setdiff(names(g), "group_id"), drop = FALSE]
    p <- p[, setdiff(names(p), "group_id"), drop = FALSE]

    # Rename probability columns back to 'prob'
    names(g)[names(g) == "prob_g"] <- "prob"
    names(p)[names(p) == "prob_p"] <- "prob"

    list(
        g_weights = g[, c(adjust_vars, time_name, "prob"), drop = FALSE],
        p_weights = p[, c(adjust_vars, "prob"), drop = FALSE]
    )
}


#' Estimate the marginalizing distribution in the observed or matched analysis-eligible
#' populations.
#'
#' @description
#' `get_observed_gp` returns marginalizing distributions  based on the original observed data
#' `get_matching_gp` returns marginalizing distributions based on the matched-analysis data
#'
#' @inheritParams get_gp
#' @return A list of two data frames containing the marginalizing distributions returned from calls to [get_gp()]
#' @keywords internal
#'
get_observed_gp <- function(data, outcome_name, trt_name,
                                            time_name, adjust_vars, tau){

    get_gp(data, outcome_name, trt_name, time_name, adjust_vars, tau)

}



#' Get empirical probability distributions g(d|x) and p(x) based on input data
#'
#' @inheritParams  nomatchVE
#' @param df A data frame representing the target population for the estimated distributions
#' * typically, population involves the "analysis-eligible" subset of a given data source
#' @param tau Numeric. Only vaccinated who are at-risk tau days after vaccination are included.
#'
#' @return A list of the following:
#' \describe{
#' \item{g_weights}{A data frame containing the covariate-conditional day probabilities.
#' Columns include `group_id` identifying the covariate group, `<time_name> specifying the day,
#' `prob` specifying the covariate-conditional day probability, and each variables in `adjust_vars`}
#' \item{p_weights}{A data frame containing the covariate probabilities.
#' Columns include `group_id` identifying the covariate group,`prob` specifying the covariate probability, and
#' each variable in `adjust_vars`}
#' }
#'
#'@keywords internal
get_gp <- function(df, outcome_name, trt_name, time_name, adjust_vars, tau){
    # Single subset operation
    gp_data <- df[(df[[outcome_name]] - df[[time_name]]) > tau & df[[trt_name]] == 1, ]

    if (nrow(gp_data) == 0L) return(list(g_weights = data.frame(), p_weights = data.frame()))

    # Create group id and merge back into g
    group_lookup <- unique(gp_data[adjust_vars])
    group_lookup <- group_lookup[do.call(order, group_lookup[adjust_vars]),]
    group_lookup$group_id <- seq_len(nrow(group_lookup))
    gp_data <- merge(gp_data, group_lookup, by = adjust_vars, all.x = TRUE)

    # Compute g_weights: P(time | covariates)
    xt <- with(gp_data, xtabs(~ group_id + gp_data[[time_name]]))
    g_prop <- prop.table(xt, 1)
    g_weights <- as.data.frame(as.table(g_prop), stringsAsFactors = FALSE)
    names(g_weights) <- c("group_id", time_name, "prob_g")
    g_weights$group_id <- as.numeric(g_weights$group_id)
    g_weights[[time_name]] <- as.numeric(g_weights[[time_name]])
    g_weights <- g_weights[g_weights$prob_g > 0, ]
    #order columns
    g_weights <- g_weights[order(g_weights$group_id, g_weights[[time_name]]),]

    # Compute p_weights: P(covariates)
    p_weights <- as.data.frame(prop.table(table(gp_data$group_id)), stringsAsFactors = FALSE)
    names(p_weights) <- c("group_id", "prob_p")
    #order columns
    p_weights$group_id <- as.numeric(p_weights$group_id)
    p_weights <- p_weights[order(p_weights$group_id),]

    # Create key for covariate values (only unique combinations)
    key <- unique(gp_data[, c("group_id", adjust_vars)])

    # Single merge for each output
    g_dist_clean <- merge(g_weights, key, by = "group_id")
    p_dist_clean <- merge(p_weights, key, by = "group_id")

    list(g_weights = g_dist_clean[, c("group_id", adjust_vars, time_name, "prob_g")],
         p_weights = p_dist_clean[, c("group_id", adjust_vars, "prob_p")])
}
