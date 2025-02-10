#' Perform n_boot bootstrap replicates
#'
#' This function is a wrapper around [one_boot_ve()] to get repeated bootstrap
#' estimates.
#'
#' @inheritParams estimate_ci
#'
#' @return {A matrix containing estimates from `n_boot` bootstrap replications. Rows
#'  represent bootstrap iterations, columns the term estimated.}
#' @export
#'

run_bootstrap <- function(data,
                          outcome_name, event_name, trt_name, time_name,
                          adjust_vars, marginalizing_dist,
                          t0, censor_time, tau,
                          boot_formula_0, boot_formula_1,
                          matched_data,
                          gp_list,
                          limit_type,
                          n_boot){

    cat("Bootstrapping...\n")
    start <- Sys.time()
    # Bootstrap using replicate function(): concise but hard to track progress
    # boot_samples <- replicate(n_boot,
    #                       one_boot_ve(data, outcome_name, event_name, trt_name, time_name, id_name,
    #                                   adjust_vars, t0, tau, formula_0, formula_1, marginalizing_dist,
    #                                   gp_list, matched_data, limit_type))
    # boot_samples <- t(boot_samples)

    # Bootstrap using for loop
    boot_samples <- matrix(NA, nrow = n_boot, ncol = 3)
    colnames(boot_samples) <- c("psi_bar_0", "psi_bar_1", "ve")
    for(i in 1:n_boot){
        #set.seed(i)
        if(i %% 50 == 0) print(i)
        boot_samples[i,] <- one_boot_ve(data = data,
                                        outcome_name = outcome_name,
                                        event_name = event_name,
                                        trt_name = trt_name,
                                        time_name = time_name,
                                        adjust_vars = adjust_vars ,
                                        marginalizing_dist = marginalizing_dist,
                                        t0 = t0,
                                        censor_time = censor_time,
                                        tau = tau,
                                        boot_formula_0 = boot_formula_0,
                                        boot_formula_1 = boot_formula_1,
                                        matched_data = matched_data,
                                        gp_list = gp_list,
                                        limit_type = limit_type)
    }

    #Boostrap run time
    end <- Sys.time()
    print(end - start)

    return(boot_samples)

}


#' Compute one bootstrap replicate of VE point estimate
#'
#' @description This function creates a bootstrapped sample and computes the
#' corresponding point estimate
#'
#' @inheritParams estimate_ci
#'
#' @return  A named numeric vector containing point estimates for
#'  `psi_bar_0, psi_bar_1, ve` in the bootstrapped sample
#'
#' @export
#'
one_boot_ve <- function(data,
                        outcome_name, event_name, trt_name, time_name,
                        adjust_vars, marginalizing_dist,
                        t0, censor_time, tau,
                        boot_formula_0,
                        boot_formula_1,
                        matched_data,
                        gp_list,
                        limit_type){

    # --------------------------------------------------------------------------
    # 1. Create bootstrapped sample(s)
    # --------------------------------------------------------------------------
    #bootstrap original data
    boot_inds <- sample(1:nrow(data), replace = TRUE)
    boot_data <-  data[boot_inds,]


    #bootstrap matched adata if needed
    if(!is.list(marginalizing_dist) && marginalizing_dist == "matched" && limit_type == "limit"){
        #print("Bootstrap matched")
        boot_match_ids <- sample(unique(matched_data$pair_id), replace = TRUE)
        boot_match_inds <- as.vector(sapply(boot_match_ids, \(pair_id) which(matched_data$pair_id == pair_id)))
        boot_match_data <- matched_data[boot_match_inds,]
        boot_match_data$pair_id <- rep(1:length(boot_match_ids), each = 2)
    }else{
        boot_match_data <- NULL
    }

    # --------------------------------------------------------------------------
    # 2. Set marginalizing distribution based on limit type
    # --------------------------------------------------------------------------
    if(is.list(marginalizing_dist)){
        limit_type <- "fixed"
        boot_marginalizing_dist <- marginalizing_dist
    }else if(limit_type == "fixed"){
        boot_marginalizing_dist <- gp_list
    }else if(limit_type == "limit"){
        boot_marginalizing_dist <- marginalizing_dist
    }else{
        stop("not a valid limit_type")
    }

    # --------------------------------------------------------------------------
    # 3. Compute VE for bootstrapped data
    # --------------------------------------------------------------------------

    boot_ve <- get_one_ve(data = boot_data ,
                          outcome_name = outcome_name,
                          event_name = event_name,
                          trt_name = trt_name,
                          time_name = time_name,
                          adjust_vars = adjust_vars ,
                          marginalizing_dist = boot_marginalizing_dist,
                          matched_dist_options = matched_dist(matched_data = boot_match_data),
                          t0 = t0,
                          censor_time = censor_time,
                          tau = tau,
                          formula_0 = boot_formula_0,
                          formula_1 = boot_formula_1,
                          return_models = FALSE,
                          return_gp_list = FALSE,
                          return_matching = FALSE)

    boot_ve$estimates

}



