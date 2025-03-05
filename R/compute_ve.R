#' Compute proposed causal VE estimand given the fitted survival models at various
#' time points of interest
#'
#' @description
#' `compute_ve()` returns the VE estimates at multiple time points. This is function is a
#' wrapper around `compute_ve_t0()`.
#'
#' @inheritParams compute_psi_d
#' @inheritParams obsve
#'
#' @return A matrix of estimates where the columns are the terms `risk_0, risk_1`, and `ve`
#' and the rows are the time points of interest.
#' @export
#'
compute_ve <- function(fit_0, fit_1, time_name, times, tau, gp_list){
    ve <- sapply(times, \(t){
        compute_ve_t0(fit_0, fit_1, time_name, t0 = t, tau, gp_list)
    })
    colnames(ve) <- times
    t(ve)
}



#' Compute proposed causal VE estimand given the fitted survival models at a single
#' time points of interest
#'
#' @inheritParams compute_psi_d
#' @inheritParams obsve
#'
#' @return A named numeric vector containing point estimates for
#'  `risk_0, risk_1`, and `ve`
#'
compute_ve_t0 <- function(fit_0, fit_1, time_name, t0, tau, gp_list){

    psi_d <- compute_psi_d(fit_0, fit_1, time_name, t0, tau, gp_list)
    psi_bar <- compute_psi_bar(psi_d, gp_list)

    # if(any(is.na(psi_bar))){
    #     cat("psi_bar: ", psi_bar, "\n")
    #     print(summary(fit_0))
    #     print(summary(fit_1))
    # }

    add_ve(psi_bar)
}

#' Calculate VE based on risk_0 and risk_1
#'
#' This function computes VE as the proportion reduction in risk due to
#' vaccination (ie. VE = 1 - risk_0 / risk_1)
#'
#' @param psi_bar {A named numeric vector containing point estimates for
#'  `risk_0, risk_1`}
#'
#' @return {A named numeric vector containing point estimates for
#'  `risk_0, risk_1`, and `ve`}
#'
#' @noRd
add_ve <- function(psi_bar){

    ve <- unname(1 - psi_bar["risk_1"]/psi_bar["risk_0"])
    c(psi_bar, "ve" = ve)
}
