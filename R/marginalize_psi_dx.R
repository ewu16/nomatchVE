#' Compute overall cumulative incidence psi_bar via marginalization
#'
#' @description This function marginalizes psi_v(t0; d, x) over estimated or pre-specified
#' marginalizing distributions
#' @param psi_d A data frame containing estimates of psi_v(t0; d, x). Data frame should
#' contain the columns `risk_0, risk_1` representing psi_0(t0; d, x) and psi_1(t0; d, x) respectively.
#' @param gp_list A list containing the marginalizing distributions
#' #TODO: assumptions on gp_list
#'
#' @return {A named numeric vector containing point estimates for
#'  `risk_0, risk_1`}
#' @export
#'
compute_psi_bar <- function(psi_d, gp_list){
   #marginalize over g*(d | x)
   #check if there are any missing covariate groups in psi_d
   psi_d$dummy <- 1
   psi_d <- merge(gp_list$g_dist, psi_d, all.x = TRUE, all.y = FALSE)

   stopifnot("psi_d missing predictions needed for marginalization" = !any(is.na(psi_d$dummy)))
   # na_risk <- which(is.na(psi_d$risk_0) | is.na(psi_d$risk_1))
   # if(length(na_risk) > 0){
   #    print(psi_d[na_risk,])
   #    return(c("risk_0" = NA, "risk_1" = NA))
   # }

   #calculate weighted mean
   psi_x <- stats::aggregate(cbind(weighted_risk_0 = psi_d$risk_0*psi_d$prob ,
                            weighted_risk_1 = psi_d$risk_1*psi_d$prob ) ~ group_name,
             data = psi_d, FUN = sum, na.action = stats::na.pass )

   #marginalize over P*(x)
   psi_x <- merge(psi_x, gp_list$p_dist)

   risk_0 <- stats::weighted.mean(psi_x$weighted_risk_0, psi_x$prob)
   risk_1 <- stats::weighted.mean(psi_x$weighted_risk_1, psi_x$prob)

   c("risk_0" = risk_0, "risk_1" = risk_1)
}


