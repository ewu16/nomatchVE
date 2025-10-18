#' Marginalize conditional cumulative incidences
#'
#' @description
#' Averages day- and covariate-specific cumulative incidences from [compute_psi_dx_t0()]
#' to produce overall cumulative incidences under each exposure. By default, performs
#' a simple average. If weights are provided for in `gp_list`, performs a weighted
#' average.
#'
#' @param psi_dx A data frame with conditional risk predictions containing:
#'   - `<adjust_vars>`: all covariate columns
#'   - `<time_name>`: vaccination day column
#'   - `psi_0_dx`: cumulative incidence under no vaccine,\eqn{\psi_0(t_0; d, x)}
#'   - `psi_1_dx`: cumulative incidence  under vaccine,\eqn{\psi_1(t_0; d, x)}
#'
#'   Typically output from [compute_psi_dx_t0()].
#'
#' @param gp_list A list of marginalizing weights with two components:
#' \describe{
#' \item{g_weights}{Data frame of covariate-conditional exposure-day probabilities \eqn{g(d \mid x)}.
#'    Must include:
#'    - `group_id`:covariate group identifier
#'    - `<time_name>`: exposure time variable
#'    - `prob_g`: probability of exposure time given the covariates
#'     - all variables in `adjust_vars`
#'  }
#' \item{p_weights}{Data frame of covariate probabilities \eqn{p(x)}.
#'     Must include:
#'     - `group_id`: covariate group identifier
#'     - `prob_p`: marginal probability of each covariate group
#'     - all variables in `adjust_vars`
#'  }
#' }
#'
#' Default is `NULL` in which case each row of `psi_dx` gets equal weight.
#'
#' @param show_intermediates Logical that only applies when `gp_list` is not `NULL;
#'  when `FALSE` (default), the function performs marginalization
#'   in a single step and returns only the overall cumulative incidences
#'   \eqn{\bar{\psi}_0(t_0)} and \eqn{\bar{\psi}_1(t_0)}.
#'   When `TRUE``, performs the two-step marginalization
#'   (first over days \eqn{d}, then covariate groups \eqn{x})
#'   and returns intermediate group-level results \eqn{\psi_v(t_0; x)}.
#'
#' @return When `show_intermediates = FALSE`, a named numeric vector with two elements:
#'   - `cuminc_0`: Overall cumulative incidence under no exposure
#'   - `cuminc_1`: Overall cumulative incidence under exposure
#'
#'   When `show_intermediates = TRUE`, a list containing both overall results (as a vector)
#'   and group-specific results (as a data frame).
#
#'
#' @seealso
#' - [compute_psi_dx_t0()] for generating `psi_dx` input
#' - [canonicalize_weights()] for creating `gp_list`
#'
#' @export
#'
marginalize_psi_dx_t0 <- function(psi_dx, gp_list = NULL, show_intermediates = FALSE){

   if(is.null(gp_list)){
      #simple marginalization- just average over predictions
      out <- c("cuminc_0" = mean(psi_dx$psi_0_dx),
               "cuminc_1" = mean(psi_dx$psi_1_dx))
      return(out)
   }

   g_weights <- gp_list$g_weights
   p_weights <- gp_list$p_weights

   if(show_intermediates == FALSE){
      # One-step marginalization - simpler and faster
      # first compute combined weight for (d,x)
      gp <- merge(gp_list$g_weights, gp_list$p_weights, suffixes = c("_g", "_p"))
      gp$prob_gp <- gp$prob_g*gp$prob_p

      #compute weighted average
      dx <- merge(psi_dx, gp)
      out <- c("cuminc_0" = sum(dx$psi_0_dx*dx$prob_gp),
               "cuminc_1" = sum(dx$psi_1_dx*dx$prob_gp))
      return(out)
   }

   # Two-step marginalization - only runs if show_intermediates = TRUE

   # Step 1: Marginalize over days within each covariate group
   dx <- merge(g_weights, psi_dx, all.x = TRUE, all.y = FALSE)
   if (anyNA(dx$psi_0_dx) | anyNA(dx$psi_1_dx)) {
      stop("psi_dx missing predictions needed for marginalization")
   }

   x <- rowsum(
      dx[c("psi_0_dx", "psi_1_dx")] * dx$prob_g,
      dx$group_id
   )
   names(x) <- c("psi_0_x", "psi_1_x")
   x$group_id <- rownames(x)
   rownames(x) <- NULL

   # Step 2: Marginalize over covariate groups
   x <- merge(p_weights, x)

   overall <- c("cuminc_0" = sum(x$psi_0_x *x$prob_p),
                "cuminc_1" = sum(x$psi_1_x *x$prob_p))

   list(psi_bar = overall, psi_x = x)

}


# marginalize_psi_dx_t0 <- function(psi_dx, gp_list){
#    #check if there are any missing covariate groups in psi_dx
#    # missing_groups <- setdiff(
#    #    unique(gp_list$g_weights$group_id),
#    #    unique(psi_dx$group_id)
#    # )
#    # if (length(missing_groups) > 0) {
#    #    stop("psi_dx missing predictions for groups: ",
#    #         paste(missing_groups, collapse = ", "))
#    # }
#    psi_dx$dummy <- 1
#    psi_dx <- merge(gp_list$g_weights, psi_dx, all.x = TRUE, all.y = FALSE)
#    stopifnot("psi_dx missing predictions needed for marginalization" = !any(is.na(psi_dx$dummy)))
#
#
#    # Step 1: Marginalize over days within each covariate group g(d|x)
#    # results in psi_0(x), psi_1(x) for each covariage group x
#    psi_x <- stats::aggregate(
#       cbind(
#          weighted_cuminc_0 = psi_dx$psi_0_dx * psi_dx$prob_g,
#          weighted_cuminc_1 = psi_dx$psi_1_dx * psi_dx$prob_g
#       ) ~ group_id,
#       data = psi_dx,
#       FUN = sum,
#       na.action = stats::na.pass
#    )
#
#    # Step 2: Marginalize over covariate groups p(x)
#    psi_x <- merge(psi_x, gp_list$p_weights)
#    overall_cuminc_0 <- stats::weighted.mean(psi_x$weighted_cuminc_0, psi_x$prob_p)
#    overall_cuminc_1 <- stats::weighted.mean(psi_x$weighted_cuminc_1, psi_x$prob_p)
#
#    c("cuminc_0" = overall_cuminc_0, "cuminc_1" = overall_cuminc_1)
# }

