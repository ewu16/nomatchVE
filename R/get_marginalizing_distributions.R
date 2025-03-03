#' Estimate the marginalizing distribution in the observed or matched analysis-eligible
#' populations.
#'
#' @description
#' `get_observed_gp` returns marginalizing distributions  based on the original observed data
#' `get_matching_gp` returns marginalizing distributions based on the matched-analysis data
#'
#' @inheritParams obsve
#' @return A list of two data frames containing the marginalizing distributions returned from calls to [get_gp()]
#' @export
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
#' @inheritParams  obsve
#' @param df A data frame representing the target population for the estimated distributions
#' * typically, population involves the "analysis-eligible" subset of a given data source
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
#'
#' @export
#'
get_gp <- function(df, outcome_name, trt_name, time_name, adjust_vars, tau){
    gp_data <-subset(df, (df[[outcome_name]] - df[[time_name]]) > tau  &
                         df[[trt_name]] == 1)

    #Define and create covariate groups
    covariate_groups <- split(gp_data, gp_data[,adjust_vars], drop = TRUE)
    key <- unique(gp_data[,adjust_vars])
    key$group_name  <- apply(key, 1, \(x) paste(trimws(x), collapse = "."))

    #Compute probabilities
    get_empirical_prob <- function(x){
        data.frame(value = sort(unique(x)),
                   prob = as.numeric(prop.table(table(x))))
    }

    g_dist <- lapply(seq_along(covariate_groups), \(i){
        group <- covariate_groups[[i]]
        prob <- cbind(names(covariate_groups)[i], get_empirical_prob(group[[time_name]]))
        names(prob) <- c("group_name", time_name, "prob")
        prob}
    ) |> stats::setNames(names(covariate_groups))

    p_dist <- sapply(covariate_groups, \(group) nrow(group)/nrow(gp_data))

    #Tidy probabilities to data frames
    g_dist_df <- do.call(rbind, c(g_dist, make.row.names = FALSE))
    p_dist_df <- data.frame(group_name = names(p_dist), prob = unname(p_dist))

    g_dist_clean <- merge(g_dist_df, key)
    p_dist_clean <- merge(p_dist_df, key)
    list(g_dist = g_dist_clean, p_dist = p_dist_clean)
}
