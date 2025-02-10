#' Simulated dataset
#'
#' Minimal data set mimicking data from an observational vaccine study.
#' Used to apply/test methods.
#'
#' @format ## `simdata`
#' A data frame with 10,000 rows and 7 columns:
#' \describe{
#'   \item{ID}{Subject identifier}
#'   \item{x1,x2}{Covariates to adjust for}
#'   \item{V}{Indicator of ever vaccinated while under follow-up}
#'   \item{D_obs}{Day of observed vaccination (relative to d0)}
#'   \item{Y}{Time to event in days (relative to d0)}
#'   \item{event}{Event indicator, event = 0 if censored}
#' }
#' @source Created from code
"simdata"
