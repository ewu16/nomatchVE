simdata_raw  <- readRDS("../../sim/sim_02/setting2/setting2_simdata.rds")
set.seed(1234)
n <- 10000
inds <- sample(1:nrow(simdata_raw), size = n, replace = FALSE)
vars <- c("ID", "x1", "x2",  "Y", "event", "V", "D_obs")

simdata <- simdata_raw[inds,vars]
simdata$ID <- 1:n
rownames(simdata) <- NULL

usethis::use_data(simdata, overwrite = TRUE)
