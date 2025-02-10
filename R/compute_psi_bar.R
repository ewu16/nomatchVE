#' Compute overall cumulative incidence psi_bar via marginalization
#'
#' @description This function marginalizes psi_v(d,x) over estimated or pre-specified
#' marginalizing distributions
#' @param psi_d A data frame containing estimates of psi_v(d,x). Data frame should
#' contain the columns `risk_0, risk_1,
#' @param gp_list A list containing the marginalizing distributions
#'
#' @return {A named numeric vector containing point estimates for
#'  `psi_bar_0, psi_bar_1`}
#' @export
#'
compute_psi_bar <- function(psi_d, gp_list){

   #marginalize over g*(d | x)
   #calculate weighted mean
   psi_x <- stats::aggregate(cbind(weighted_risk_0 = psi_d$risk_0*psi_d$prob ,
                            weighted_risk_1 = psi_d$risk_1*psi_d$prob ) ~ group_name,
             data = psi_d, FUN = sum, na.action = stats::na.pass )

   #marginalize over P*(x)
   psi_x <- merge(psi_x, gp_list$p_dist)

   psi_bar_0 <- stats::weighted.mean(psi_x$weighted_risk_0, psi_x$prob)
   psi_bar_1 <- stats::weighted.mean(psi_x$weighted_risk_1, psi_x$prob)

   c("psi_bar_0" = psi_bar_0, "psi_bar_1" = psi_bar_1)
}


#' Calculate VE based on psi_bar_0 and psi_bar_1
#'
#' This function computes VE as the proportion reduction in risk due to
#' vaccination (ie. VE = 1 - psi_bar_0 / psi_bar_1)
#'
#' @param psi_bar {A named numeric vector containing point estimates for
#'  `psi_bar_0, psi_bar_1`}
#'
#' @return {A named numeric vector containing point estimates for
#'  `psi_bar_0, psi_bar_1`, and `ve`}
#' @export
#'
add_ve <- function(psi_bar){

   ve <- unname(1 - psi_bar["psi_bar_1"]/psi_bar["psi_bar_0"])
   c(psi_bar, "ve" = ve)
}


