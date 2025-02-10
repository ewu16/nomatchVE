
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`obsve`

> Estimate vaccine effectiveness in observational studies: a matching
> alternative

<!-- badges: start -->

<!-- badges: end -->

## Description

The `obsve` package uses a G-computation style estimator to compute
vaccine efficacy from observational vaccine studies. The proposed
estimator tends to produce similar point estimates as matching-based
estimators but is more efficient.

## Installation

You can install the development version of `obsve` like so:

``` r
# TODO: not yet available on Github 
# install.packages("devtools")
devtools::install_github("ewu16/obsve")
```

## Example

This minimal example shows how to use `obsve` to obtain cumulative
incidence and vaccine effectiveness estimates in a simple simulated data
set.

``` r
library(obsve)
# Set seed for reproducibility 
set.seed(1234)

# ------------------------------------------------------------------------------
# Example data
head(simdata)
#>   ID x1 x2 V D_obs   Y event
#> 1  1  1  7 1     2  92     0
#> 2  2  0  7 0    NA 210     0
#> 3  3  0 11 1    35 210     0
#> 4  4  0 10 1     6 210     0
#> 5  5  1 11 0    NA 210     0
#> 6  6  1  7 0    NA  90     0
summary(simdata)
#>        ID              x1               x2               V         
#>  Min.   :    1   Min.   :0.0000   Min.   : 5.000   Min.   :0.0000  
#>  1st Qu.: 2501   1st Qu.:0.0000   1st Qu.: 6.000   1st Qu.:0.0000  
#>  Median : 5000   Median :0.0000   Median : 8.000   Median :0.0000  
#>  Mean   : 5000   Mean   :0.4989   Mean   : 8.023   Mean   :0.4112  
#>  3rd Qu.: 7500   3rd Qu.:1.0000   3rd Qu.:10.000   3rd Qu.:1.0000  
#>  Max.   :10000   Max.   :1.0000   Max.   :11.000   Max.   :1.0000  
#>                                                                    
#>      D_obs              Y           event       
#>  Min.   :  1.00   Min.   :  1   Min.   :0.0000  
#>  1st Qu.: 11.00   1st Qu.:174   1st Qu.:0.0000  
#>  Median : 18.00   Median :210   Median :0.0000  
#>  Mean   : 25.78   Mean   :178   Mean   :0.1007  
#>  3rd Qu.: 32.00   3rd Qu.:210   3rd Qu.:0.0000  
#>  Max.   :206.00   Max.   :210   Max.   :1.0000  
#>  NA's   :5888

# ------------------------------------------------------------------------------
# 1. Set input parameters
outcome_name <- "Y"
event_name <- "event"
trt_name <- "V"
time_name <- "D_obs"
adjust_vars <- c("x1", "x2")

t0 <- 90
censor_time <- 90
tau <- 14
ci_type <- "wald"
limit_type <- "fixed"
n_boot <- 10
alpha <- .05

# ------------------------------------------------------------------------------
# 2. Compute VE estimand at t0 

## 2a. Use observed distribution of vaccination dates/covariates

fit1 <- obsve(data = simdata,
              outcome_name = outcome_name,
              event_name = event_name,
              trt_name = trt_name,
              time_name = time_name, 
              adjust_vars = adjust_vars,
              marginalizing_dist = "observed",
              t0 = t0,
              censor_time = censor_time, 
              tau = tau,
              ci_type = ci_type,
              limit_type =  limit_type,
              n_boot = n_boot, 
              alpha = alpha)
#> Bootstrapping...
#> Time difference of 0.2809708 secs
              
             
fit1$estimates
#>             estimate wald_lower wald_upper    wald_sd     boot_sd
#> psi_bar_0 0.05857293 0.05481528 0.06257112 0.03586146 0.002395859
#> psi_bar_1 0.03520708 0.03029040 0.04088816 0.07933783 0.002863597
#> ve        0.39891903 0.29342765 0.48866052 0.08249945 0.047365506



## 2b. Use distribution of vaccination dates/covariates from matched data set 

### Create matched dataset to pass as additional argument to `obsve()`
id_name <- "ID"
matching_vars <- adjust_vars

matched_cohort <- match_rolling_cohort(data = simdata,
                                       outcome_name = outcome_name,
                                       trt_name = trt_name,
                                       time_name = time_name,
                                       id_name = id_name,
                                       matching_vars = adjust_vars,
                                       replace = FALSE,
                                       seed = 5678)

matched_data <- matched_cohort[[1]]

fit2 <- obsve(data = simdata,
              outcome_name = outcome_name,
              event_name = event_name,
              trt_name = trt_name,
              time_name = time_name,
              adjust_vars = adjust_vars,
              marginalizing_dist = "matched",
              matched_dist_options = matched_dist(matched_data = matched_data,
                                                  pair_censoring = TRUE),
              t0 = t0,
              censor_time = censor_time, 
              tau = tau,
              ci_type = ci_type,
              limit_type =  limit_type,
              n_boot = n_boot, 
              alpha = alpha)
#> Bootstrapping...
#> Time difference of 0.3211391 secs


## Either choice of marginalizing distribution leads to similar estimates
fit1$estimates
#>             estimate wald_lower wald_upper    wald_sd     boot_sd
#> psi_bar_0 0.05857293 0.05481528 0.06257112 0.03586146 0.002395859
#> psi_bar_1 0.03520708 0.03029040 0.04088816 0.07933783 0.002863597
#> ve        0.39891903 0.29342765 0.48866052 0.08249945 0.047365506

fit2$estimates
#>             estimate wald_lower wald_upper    wald_sd     boot_sd
#> psi_bar_0 0.05810865 0.05397950 0.06253277 0.03983969 0.002673693
#> psi_bar_1 0.03533592 0.02989989 0.04171778 0.08809592 0.003307776
#> ve        0.39189904 0.24734115 0.50869272 0.10881381 0.062919660

# ------------------------------------------------------------------------------
# 3. Compute VE at additional timepoints if desired 

times <- c(30, 60, 90)
results1 <- timepoints(fit1, times = times)
#> Timepoint: 30 
#> Bootstrapping...
#> Time difference of 0.2827771 secs
#> Timepoint: 60 
#> Bootstrapping...
#> Time difference of 0.3260989 secs
results1$estimates
#> $`30`
#>              estimate  wald_lower wald_upper    wald_sd      boot_sd
#> psi_bar_0 0.011638913 0.010809067 0.01253166 0.03816803 0.0004618664
#> psi_bar_1 0.006211499 0.003714973 0.01036827 0.26354558 0.0021087436
#> ve        0.466316230 0.049543685 0.70033513 0.29446395 0.1978513190
#> 
#> $`60`
#>             estimate wald_lower wald_upper    wald_sd     boot_sd
#> psi_bar_0 0.03790915 0.03426256 0.04192700 0.05353284 0.002338379
#> psi_bar_1 0.02282077 0.01895132 0.02745817 0.09681284 0.002394330
#> ve        0.39801415 0.22411067 0.53293989 0.12947991 0.075504986
#> 
#> $`90`
#>             estimate wald_lower wald_upper    wald_sd     boot_sd
#> psi_bar_0 0.05857293 0.05481528 0.06257112 0.03586146 0.002395859
#> psi_bar_1 0.03520708 0.03029040 0.04088816 0.07933783 0.002863597
#> ve        0.39891903 0.29342765 0.48866052 0.08249945 0.047365506

results2 <- timepoints(fit2, times = times)
#> Timepoint: 30 
#> Bootstrapping...
#> Time difference of 0.2718439 secs
#> Timepoint: 60 
#> Bootstrapping...
#> Time difference of 0.277036 secs


# ------------------------------------------------------------------------------
# 4. Compare results with matching estimator

fit_matching <-matching_ve(matched_data = matched_data,
                           outcome_name = outcome_name,
                           event_name = event_name,
                           trt_name = trt_name,
                           time_name = time_name,
                           adjust = adjust_vars,
                           times = times,
                           censor_time = censor_time,
                           tau = tau,
                           pair_censoring = TRUE,
                           separate = TRUE,
                           ci_type = ci_type,
                           limit_type = limit_type,
                           data = simdata,
                           id_name = id_name,
                           matching_vars = matching_vars,
                           n_boot = n_boot,
                           alpha = alpha) 
#> Bootstrapping...
#> Time difference of 1.067213 secs
   

## Proposed and matching based estimators have similar point estimates
##  Proposed has narrower confidence intervals 
results1$estimates
#> $`30`
#>              estimate  wald_lower wald_upper    wald_sd      boot_sd
#> psi_bar_0 0.011638913 0.010809067 0.01253166 0.03816803 0.0004618664
#> psi_bar_1 0.006211499 0.003714973 0.01036827 0.26354558 0.0021087436
#> ve        0.466316230 0.049543685 0.70033513 0.29446395 0.1978513190
#> 
#> $`60`
#>             estimate wald_lower wald_upper    wald_sd     boot_sd
#> psi_bar_0 0.03790915 0.03426256 0.04192700 0.05353284 0.002338379
#> psi_bar_1 0.02282077 0.01895132 0.02745817 0.09681284 0.002394330
#> ve        0.39801415 0.22411067 0.53293989 0.12947991 0.075504986
#> 
#> $`90`
#>             estimate wald_lower wald_upper    wald_sd     boot_sd
#> psi_bar_0 0.05857293 0.05481528 0.06257112 0.03586146 0.002395859
#> psi_bar_1 0.03520708 0.03029040 0.04088816 0.07933783 0.002863597
#> ve        0.39891903 0.29342765 0.48866052 0.08249945 0.047365506

fit_matching$estimates
#> $`30`
#>           estimate  wald_lower  wald_upper    wald_sd           sd
#> risk_0 0.010864854 0.008219022 0.014350098 0.14375433 0.0017218990
#> risk_1 0.005403252 0.004570544 0.006386699 0.08582084 0.0004586441
#> ve     0.502685234 0.362570138 0.612001271 0.12664574 0.0558723196
#> 
#> $`60`
#>          estimate wald_lower wald_upper   wald_sd          sd
#> risk_0 0.04053786 0.03337100 0.04916561 0.1030585 0.004661797
#> risk_1 0.02094066 0.01668776 0.02624850 0.1180389 0.002774174
#> ve     0.48342944 0.27644009 0.63120518 0.1719274 0.096222447
#> 
#> $`90`
#>          estimate wald_lower wald_upper    wald_sd          sd
#> risk_0 0.05650257 0.04895697 0.06513152 0.07720051 0.004942439
#> risk_1 0.03451473 0.02700801 0.04401351 0.12908277 0.005115813
#> ve     0.38914763 0.12566488 0.57322930 0.18296685 0.118559140
```
