
# library(BayesianTools)
library(cmdstanr)

results = readRDS("./toyPara_GPUs_nonInformative_noProcessError_23-08-2021_12h22.rds")

mcmc_pairs(results$draws(), regex_pars = c("Error", "intercept"))

