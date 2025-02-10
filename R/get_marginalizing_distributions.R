#get_marginalizing_distributions

#' Estimate the marginalizing distribution in the observed or matched analysis-eligible
#' populations.
#'
#' @inheritParams obsve
#' @param matched_adata Needed when marginalizing_dist_type = "matched". This is
#' a data frame representing the matched analysis-eligible population
#' @param marginalizing_dist_type Character string describing the type of estimated
#'   day/covariate distributions to marginalize over: `"observed"` or "`matched`"
#' @return A list of two data frames containing the marginalizing distributions
#' @export
#'
get_marginalizing_distributions <- function(data, outcome_name, trt_name,
                                            time_name, adjust_vars, tau,
                                            marginalizing_dist_type,
                                            matched_adata = NULL){

    stopifnot("not a valid marginalizing_dist_type" =
                  marginalizing_dist_type %in% c("observed", "matched"))


    if(marginalizing_dist_type == "observed"){
        gp_list <-  get_gp(data, outcome_name, trt_name, time_name, adjust_vars, tau)

    }else if(marginalizing_dist_type == "matched" & !is.null(matched_adata)){
        stopifnot("must provide matched adata" = !is.null(matched_adata))

        gp_list <- get_gp(df = matched_adata,
                          outcome_name = paste0(outcome_name, "_matched"),
                          trt_name = paste0(trt_name, "_d"),
                          time_name = "d",
                          adjust_vars = adjust_vars,
                          tau = tau)
        #rename d variable to original time_name
        names(gp_list$g_dist)[which(names(gp_list$g_dist ) == "d")] <- time_name
    }

    return(gp_list)


}


#' Get empirical probability distributions g(d|x) and P(x) based on input data
#'
#' @inheritParams  obsve
#' @param df A data frame representing the target population for the estimated distributions
#' * typically, population involves the "analysis-eligible" subset of a given data source
#'
#' @return A list of two data frames containing the estimated probabilities
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
