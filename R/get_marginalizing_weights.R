#' Estimate observed marginalizing distributions for G-computation
#'
#' @description
#' Computes empirical probability distributions used to marginalize predictions over
#' observed exposure timing and covariate patternss:
#' - \eqn{g(d \mid x)}: the distribution of exposure times within each covariate subgroup
#'  of exposed individuals who remain at risk `tau` days after exposure.
#' - \eqn{p(x)}: the distribution of covariates among exposed individuals
#'   who remain at risk `tau` days after exposure.
#'
#' Most users do not need to call this function directly — it is called
#' internally by [nomatchVE()]. It can be useful for:
#' - Creating or understanding custom marginalizing distributions
#' - Diagnosing weighting or estimation issues
#'
#' @inheritParams get_gp
#'
#' @return A list with two data frames:
#' - `g_dist`: covariate-conditional exposure time  probabilities (\eqn{g(d \mid x)})
#' - `p_dist`: covariate probabilities (\eqn{p(x)})
#'
#' Each data frame includes a `group_name` column identifying the covariate
#' subgroup, the columns for variables in `adjust_vars`, and a `prob` column
#' with empirical probabilities. `g_dist` additionally includes a `<time_name>`
#' column for exposure timing.
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

    get_gp(data, outcome_name, trt_name, time_name, adjust_vars, tau)

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

#' @rdname get_observed_gp
#' @param matched_adata For `get_matching_gp`, a data frame representing the matched analysis-eligible population
get_matching_gp <- function(matched_adata, outcome_name, trt_name,
                                         time_name, adjust_vars, tau){
    gp_list <- get_gp(df = matched_adata,
                      outcome_name = paste0(outcome_name, "_matched"),
                      trt_name = paste0(trt_name, "_d"),
                      time_name = "d",
                      adjust_vars = adjust_vars,
                      tau = tau)
    #rename d variable to original time_name
    names(gp_list$g_dist)[which(names(gp_list$g_dist ) == "d")] <- time_name
    gp_list
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
#' \item{g_dist}{A data frame containing the covariate-conditional day probabilities.
#' Columns include `group_name` identifying the covariate group, `<time_name> specifying the day,
#' `prob` specifying the covariate-conditional day probability, and each variables in `adjust_vars`}
#' \item{p_dist}{A data frame containing the covariate probabilities.
#' Columns include `group_name` identifying the covariate group,`prob` specifying the covariate probability, and
#' each variable in `adjust_vars`}
#' }
#'
#'@keywords internal
get_gp <- function(df, outcome_name, trt_name, time_name, adjust_vars, tau){
    # Single subset operation
    gp_data <- df[(df[[outcome_name]] - df[[time_name]]) > tau & df[[trt_name]] == 1, ]

    if (nrow(gp_data) == 0L) return(list(g_dist = data.frame(), p_dist = data.frame()))

    # Create group_name- efficient because interaction() is vectorized
    gp_data$group_name <- do.call(interaction, c(gp_data[adjust_vars], sep = ".", drop = TRUE))

    # Compute g_dist: P(time | covariates)
    xt <- with(gp_data, xtabs(~ group_name + gp_data[[time_name]]))
    g_prop <- prop.table(xt, 1)
    g_dist <- as.data.frame(as.table(g_prop), stringsAsFactors = FALSE)
    names(g_dist) <- c("group_name", time_name, "prob")
    g_dist[[time_name]] <- as.numeric(g_dist[[time_name]])
    g_dist <- g_dist[g_dist$prob > 0, ]
    #order columns
    g_dist <- g_dist[order(g_dist$group_name, g_dist[[time_name]]),]

    # Compute p_dist: P(covariates)
    p_dist <- as.data.frame(prop.table(table(group_name)))
    names(p_dist) <- c("group_name", "prob")
    #order columns
    p_dist <- p_dist[order(p_dist$group_name),]

    # Create key for covariate values (only unique combinations)
    key <- unique(gp_data[, c("group_name", adjust_vars)])

    # Single merge for each output
    g_dist_clean <- merge(g_dist, key, by = "group_name")
    p_dist_clean <- merge(p_dist, key, by = "group_name")

    list(g_dist = g_dist_clean[, c("group_name", adjust_vars, time_name, "prob")],
         p_dist = p_dist_clean[, c("group_name", adjust_vars, "prob")])
}
