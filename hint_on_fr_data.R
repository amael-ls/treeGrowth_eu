
#### Aim of prog: To have a look on the French data and get some hints on the priors
# Remember that the gamma distrib can be reparameterised by shape = mean^2/var, rate = mean/var

rm(list = ls())
graphics.off()

library(data.table)
library(bayesplot)
library(cmdstanr)
library(actuar)

options(max.print = 500)

#### Load data
trees = readRDS("../Inventories/FR IFN/trees_forest.rds")

trees[, dbh_5_in_mm := dbh_in_mm + increment_5_yrs]
trees[, ratio := dbh_5_in_mm/dbh_in_mm]

hist(trees[, ratio])
hist(trees[, dbh_in_mm])

if (!file.exists("./ratio.jpg"))
{
	jpeg("ratio.jpg", quality = 100, height = 1080, width = 1080)
	plot(trees[, dbh_in_mm], trees[, ratio], pch = 19, xlab = "Diameter t0", ylab = "ratio diameter t1/diameter t0",
		col = "#32323210")
	dev.off()
}

# x_try = 500
# y_try = 1.05
# alpha = 0.2

# fct = function(x)
# {
# 	beta = 1/x_try * log((y_try - 1)/alpha)
# 	return (1 + alpha*exp(beta*x))
# }

# curve(fct(x), from = 0, to = 1500, lwd = 2, col = "red", add = TRUE)

# plot(trees[, circumference], trees[, increment_5_yrs], pch = 19) # Gives nothing clear


#### Quick stan model on the ratio
model = cmdstan_model("ratio.stan")
maxIter = 4e3
n_chains = 3

stanData = list(N = trees[speciesName_fr == "Tilleul à grandes feuilles", .N],
	Y = trees[speciesName_fr == "Tilleul à grandes feuilles", ratio],
	dbh = trees[speciesName_fr == "Tilleul à grandes feuilles", dbh_in_mm])

results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 200, chains = n_chains,
	iter_warmup = maxIter/2, iter_sampling = maxIter/2, save_warmup = TRUE,
	max_treedepth = 11, adapt_delta = 0.95)

results$cmdstan_diagnose()

## Alpha term
plot_title = ggplot2::ggtitle("Posterior distribution alpha", "with medians and 80% intervals")
mcmc_areas(results$draws("alpha"), prob = 0.8) + plot_title

plot_title = ggplot2::ggtitle("Traces for alpha")
mcmc_trace(results$draws("alpha")) + plot_title

## Beta term
plot_title = ggplot2::ggtitle("Posterior distribution beta", "with medians and 80% intervals")
mcmc_areas(results$draws("beta"), prob = 0.8) + plot_title

plot_title = ggplot2::ggtitle("Traces for beta")
mcmc_trace(results$draws("beta")) + plot_title

## Epsilon term
plot_title = ggplot2::ggtitle("Posterior distribution epsilon", "with medians and 80% intervals")
mcmc_areas(results$draws("epsilon"), prob = 0.8) + plot_title

plot_title = ggplot2::ggtitle("Traces for epsilon")
mcmc_trace(results$draws("epsilon")) + plot_title

#### Check-up plots with pointwise values (mean). _ptw for pointwise
dbh = trees[speciesName_fr == "Tilleul à grandes feuilles", dbh_in_mm]

# results$print(variables = c("alpha", "beta", "epsilon"), digits = 10)

alpha_ptw = mean(results$draws(variables = c("alpha"), format = "draws_df")$alpha)
beta_ptw = mean(results$draws(variables = c("beta"), format = "draws_df")$beta)
gamma_ptw = mean(results$draws(variables = c("gamma"), format = "draws_df")$gamma)
epsilon_ptw = mean(results$draws(variables = c("epsilon"), format = "draws_df")$epsilon)

fct = function(x)
	return (1 + alpha_ptw*exp(-beta_ptw * x))

sigma = epsilon_ptw/dbh^gamma_ptw
n_rep = 100

dt = data.table(rep_dbh = rep(dbh, each = n_rep), sampled = numeric(n_rep * length(dbh)))

for (i in 1:length(dbh))
	dt[((i - 1)*n_rep + 1):(i*n_rep), sampled := rgamma(n_rep, shape = fct(dbh[i])^2/sigma[i], rate = fct(dbh[i])/sigma[i])]

plot(dbh, trees[speciesName_fr == "Tilleul à grandes feuilles", ratio], ylim = c(0.75, 1.5),
	pch = 19, xlab = "Diameter t0", ylab = "ratio diameter t1/diameter t0", col = "#32323222")
curve(fct(x), add = TRUE, col = "#CD212A", from = 0, to = 1500, lwd = 2)
points(dt[, rep_dbh], dt[, sampled], pch = 19, col = "#4DA6FF03")

index_min = which.min(dbh)
sampled = rgamma(1e4, shape = fct(dbh[index_min])^2/sigma[index_min], rate = fct(dbh[index_min])/sigma[index_min])


#?####################################################################
#*###########    MODEL USING A LOGNORMAL DISTRIBUTION    #############
#?####################################################################
#### Quick stan model on the ratio
model_logN = cmdstan_model("ratio_logN.stan")
maxIter = 4e3
n_chains = 3

stanData = list(N = trees[speciesName_fr == "Tilleul à grandes feuilles", .N],
	Y = trees[speciesName_fr == "Tilleul à grandes feuilles", ratio],
	dbh = trees[speciesName_fr == "Tilleul à grandes feuilles", dbh_in_mm])

results_logN = model_logN$sample(data = stanData, parallel_chains = n_chains, refresh = 200, chains = n_chains,
	iter_warmup = maxIter/2, iter_sampling = maxIter/2, save_warmup = TRUE,
	max_treedepth = 11, adapt_delta = 0.95)

results_logN$cmdstan_diagnose()

## Alpha term
plot_title = ggplot2::ggtitle("Posterior distribution alpha", "with medians and 80% intervals")
mcmc_areas(results_logN$draws("alpha"), prob = 0.8) + plot_title

plot_title = ggplot2::ggtitle("Traces for alpha")
mcmc_trace(results_logN$draws("alpha")) + plot_title

## Beta term
plot_title = ggplot2::ggtitle("Posterior distribution beta", "with medians and 80% intervals")
mcmc_areas(results_logN$draws("beta"), prob = 0.8) + plot_title

plot_title = ggplot2::ggtitle("Traces for beta")
mcmc_trace(results_logN$draws("beta")) + plot_title

## Epsilon term
plot_title = ggplot2::ggtitle("Posterior distribution epsilon", "with medians and 80% intervals")
mcmc_areas(results_logN$draws("epsilon"), prob = 0.8) + plot_title

plot_title = ggplot2::ggtitle("Traces for epsilon")
mcmc_trace(results_logN$draws("epsilon")) + plot_title

#### Check-up plots with pointwise values (mean). _ptw for pointwise
dbh = trees[speciesName_fr == "Tilleul à grandes feuilles", dbh_in_mm]

alpha_ptw = mean(results_logN$draws(variables = c("alpha"), format = "draws_df")$alpha)
beta_ptw = mean(results_logN$draws(variables = c("beta"), format = "draws_df")$beta)
epsilon_ptw = mean(results_logN$draws(variables = c("epsilon"), format = "draws_df")$epsilon)

fct = function(x)
	return (1 + alpha_ptw*exp(-beta_ptw * x))

sigma = epsilon_ptw/(dbh + 1)
n_rep = 100

dt = data.table(rep_dbh = rep(dbh, each = n_rep), sampled = numeric(n_rep * length(dbh)))

for (i in 1:length(dbh))
	dt[((i - 1)*n_rep + 1):(i*n_rep), sampled := rlnorm(n_rep, meanlog = fct(dbh[i]) - sigma[i]^2/2, sdlog = sigma[i])]

plot(dbh, trees[speciesName_fr == "Tilleul à grandes feuilles", ratio],
	pch = 19, xlab = "Diameter t0", ylab = "ratio diameter t1/diameter t0", col = "#32323210")
curve(fct(x), add = TRUE, col = "#CD212A", from = 0, to = 1500, lwd = 2)
points(dt[1:100, rep_dbh], dt[1:100, sampled], pch = 19, col = "#4DA6FF10")

index_min = which.min(dbh)
sampled = rlnorm(1e4, meanlog = fct(dbh[index_min]) - sigma[index_min]^2/2, sdlog = sigma[index_min])

mean(sampled)

plot(rep(dbh[index_min], length(sampled)), sampled)

tt = rlnorm(1e6, 0.1, 0.3)
mean(log(tt))

exp(0.1 + 0.3^2/2)
mean(tt)

sd(log(tt))
sd(tt)
(exp(0.3^2) - 1)*exp(2*0.1 + 0.3^2)



#?######################################################################
#*###########    MODEL USING A GAMMA DISTRIB + FLEX AVG    #############
#?######################################################################
#### Quick stan model on the ratio
model = cmdstan_model("ratio_flexAvg.stan")
maxIter = 4e3
n_chains = 3

stanData = list(N = trees[speciesName_fr == "Tilleul à grandes feuilles", .N],
	Y = trees[speciesName_fr == "Tilleul à grandes feuilles", ratio],
	dbh = trees[speciesName_fr == "Tilleul à grandes feuilles", dbh_in_mm])

results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 200, chains = n_chains,
	iter_warmup = maxIter/2, iter_sampling = maxIter/2, save_warmup = TRUE,
	max_treedepth = 11, adapt_delta = 0.95)

results$cmdstan_diagnose()

## Alpha term
plot_title = ggplot2::ggtitle("Posterior distribution alpha", "with medians and 80% intervals")
mcmc_areas(results$draws("alpha"), prob = 0.8) + plot_title

plot_title = ggplot2::ggtitle("Traces for alpha")
mcmc_trace(results$draws("alpha")) + plot_title

## Beta term
plot_title = ggplot2::ggtitle("Posterior distribution beta", "with medians and 80% intervals")
mcmc_areas(results$draws("beta"), prob = 0.8) + plot_title

plot_title = ggplot2::ggtitle("Traces for beta")
mcmc_trace(results$draws("beta")) + plot_title

## Epsilon term
plot_title = ggplot2::ggtitle("Posterior distribution epsilon", "with medians and 80% intervals")
mcmc_areas(results$draws("epsilon"), prob = 0.8) + plot_title

plot_title = ggplot2::ggtitle("Traces for epsilon")
mcmc_trace(results$draws("epsilon")) + plot_title

#### Check-up plots with pointwise values (mean). _ptw for pointwise
dbh = trees[speciesName_fr == "Tilleul à grandes feuilles", dbh_in_mm]

# results$print(variables = c("alpha", "beta", "epsilon"), digits = 10)

alpha_ptw = mean(results$draws(variables = c("alpha"), format = "draws_df")$alpha)
beta_ptw = mean(results$draws(variables = c("beta"), format = "draws_df")$beta)
gamma_ptw = mean(results$draws(variables = c("gamma"), format = "draws_df")$gamma)
epsilon_ptw = mean(results$draws(variables = c("epsilon"), format = "draws_df")$epsilon)

fct = function(x)
	return (1 + alpha_ptw/x^beta_ptw)
	# return (1 + alpha_ptw*exp(-beta_ptw * x))

sigma = epsilon_ptw/dbh^gamma_ptw
# sigma = epsilon_ptw/(dbh + 1)
n_rep = 100

dt = data.table(rep_dbh = rep(dbh, each = n_rep), sampled = numeric(n_rep * length(dbh)))

for (i in 1:length(dbh))
	dt[((i - 1)*n_rep + 1):(i*n_rep), sampled := rgamma(n_rep, shape = fct(dbh[i])^2/sigma[i], rate = fct(dbh[i])/sigma[i])]

plot(dbh, trees[speciesName_fr == "Tilleul à grandes feuilles", ratio], ylim = c(0.75, 1.5),
	pch = 19, xlab = "Diameter t0", ylab = "ratio diameter t1/diameter t0", col = "#323232AA")
curve(fct(x), add = TRUE, col = "#CD212A", from = 0, to = 1500, lwd = 2)
points(dt[, rep_dbh], dt[, sampled], pch = 19, col = "#4DA6FF03")

index_min = which.min(dbh)
sampled = rgamma(1e4, shape = fct(dbh[index_min])^2/sigma[index_min], rate = fct(dbh[index_min])/sigma[index_min])



#?#################################################################
#*###########    MODEL USING A PARETO DISTRIBUTION    #############
#?#################################################################
#### Quick stan model on the ratio
model = cmdstan_model("ratio_pareto.stan")
maxIter = 4e3
n_chains = 3

stanData = list(N = trees[speciesName_fr == "Tilleul à grandes feuilles", .N],
	Y = trees[speciesName_fr == "Tilleul à grandes feuilles", ratio],
	dbh = trees[speciesName_fr == "Tilleul à grandes feuilles", dbh_in_mm])

results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 200, chains = n_chains,
	iter_warmup = maxIter/2, iter_sampling = maxIter/2, save_warmup = TRUE,
	max_treedepth = 11, adapt_delta = 0.95)

results$cmdstan_diagnose()

## A term
plot_title = ggplot2::ggtitle("Posterior distribution a", "with medians and 80% intervals")
mcmc_areas(results$draws("a"), prob = 0.8) + plot_title

plot_title = ggplot2::ggtitle("Traces for a")
mcmc_trace(results$draws("a")) + plot_title

# B term
plot_title = ggplot2::ggtitle("Posterior distribution b", "with medians and 80% intervals")
mcmc_areas(results$draws("b"), prob = 0.8) + plot_title

plot_title = ggplot2::ggtitle("Traces for b")
mcmc_trace(results$draws("b")) + plot_title

## Epsilon term
plot_title = ggplot2::ggtitle("Posterior distribution epsilon", "with medians and 80% intervals")
mcmc_areas(results$draws("epsilon"), prob = 0.8) + plot_title

plot_title = ggplot2::ggtitle("Traces for epsilon")
mcmc_trace(results$draws("epsilon")) + plot_title

## Gamma term
plot_title = ggplot2::ggtitle("Posterior distribution gamma", "with medians and 80% intervals")
mcmc_areas(results$draws("gamma"), prob = 0.8) + plot_title

plot_title = ggplot2::ggtitle("Traces for gamma")
mcmc_trace(results$draws("gamma")) + plot_title

#### Check-up plots with pointwise values (mean). _ptw for pointwise
dbh = trees[speciesName_fr == "Tilleul à grandes feuilles", dbh_in_mm]

# results$print(variables = c("alpha", "beta", "epsilon"), digits = 10)

a_ptw = mean(results$draws(variables = c("a"), format = "draws_df")$a)
b_ptw = mean(results$draws(variables = c("b"), format = "draws_df")$b)
gamma_ptw = median(results$draws(variables = c("gamma"), format = "draws_df")$gamma)
epsilon_ptw = median(results$draws(variables = c("epsilon"), format = "draws_df")$epsilon)

fct = function(x)
	return (1 + a_ptw/x^b_ptw)
	# return (1 + alpha_ptw*exp(-beta_ptw * x))

n_rep = 100

dt = data.table(rep_dbh = rep(dbh, each = n_rep), sampled = numeric(n_rep * length(dbh)))

for (i in 1:length(dbh))
{
	var_i = epsilon_ptw/dbh[i]^gamma_ptw;
	avg_i = fct(dbh[i]);
	scale_i = (1 - avg_i)*(1 - 2*avg_i + avg_i^2 + var_i)/(1 - 2*avg_i + avg_i^2 - var_i); # lambda in stan
	shape_i = 2*var_i/(-1 + 2*avg_i + avg_i^2 + var_i); # alpha in stan
	dt[((i - 1)*n_rep + 1):(i*n_rep), sampled := rpareto2(n_rep, min = 1, shape = shape_i, scale = scale_i)]
}

plot(dbh, trees[speciesName_fr == "Tilleul à grandes feuilles", ratio], ylim = c(0.75, 1.5),
	pch = 19, xlab = "Diameter t0", ylab = "ratio diameter t1/diameter t0", col = "#323232AA")
curve(fct(x), add = TRUE, col = "#CD212A", from = 0, to = 1500, lwd = 2)
points(dt[, rep_dbh], dt[, sampled], pch = 19, col = "#4DA6FF03")

index_min = which.min(dbh)
var_i = epsilon_ptw/dbh[index_min]^gamma_ptw;
avg_i = fct(dbh[index_min]);
scale_i = (1 - avg_i)*(1 - 2*avg_i + avg_i^2 + var_i)/(1 - 2*avg_i + avg_i^2 - var_i); # lambda in stan
shape_i = 2*var_i/(-1 + 2*avg_i + avg_i^2 + var_i); # alpha in stan
sim = rpareto2(n_rep, min = 1, shape = shape_i, scale = scale_i)

plot(rep(dbh[index_min], n_rep), sim)
points(rep(dbh[index_min], 61), trees[speciesName_fr == "Tilleul à grandes feuilles" & dbh_in_mm == min(dbh), ratio], col = "red")

#?###################################################################################
#*###########    MODEL USING A PARETO DISTRIBUTION WITH USUAL PARAMS    #############
#?###################################################################################
#### Tool function
## Initiate Y_gen with reasonable value (by default, stan would generate them between 0 and 2---constraint Y_gen > 0)
init_fun = function(...)
{
	providedArgs = list(...)
	requiredArgs = c("dbh")
	if (!all(requiredArgs %in% names(providedArgs)))
		stop("You must provide chain_id, Yobs, and years_indiv")

	# lambda0_init = rnorm(n = 1, mean = mean(dbh), sd = 3)
	lambda1_init = sd(dbh)

	# return(list(lambda_0 = lambda0_init, lambda_1 = lambda1_init))
	return(list(lambda_1 = lambda1_init))
}

#### Quick stan model on the ratio
model = cmdstan_model("ratio_pareto_noMoment.stan")
maxIter = 4e3
n_chains = 3

# init_lambda = lapply(1:n_chains, init_fun, dbh = trees[speciesName_fr == "Tilleul à grandes feuilles", dbh_in_mm])

stanData = list(N = trees[speciesName_fr == "Tilleul à grandes feuilles", .N],
	Y = trees[speciesName_fr == "Tilleul à grandes feuilles", ratio],
	dbh = trees[speciesName_fr == "Tilleul à grandes feuilles", dbh_in_mm])

results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 200, chains = n_chains,
	iter_warmup = maxIter/2, iter_sampling = maxIter/2, save_warmup = TRUE, # init = init_lambda,
	max_treedepth = 11, adapt_delta = 0.95)

results$cmdstan_diagnose()

# Alpha_0
plot_title = ggplot2::ggtitle("Posterior distribution alpha_0", "with medians and 80% intervals")
mcmc_areas(results$draws("alpha_0"), prob = 0.8) + plot_title

plot_title = ggplot2::ggtitle("Traces for alpha_0")
mcmc_trace(results$draws("alpha_0")) + plot_title

# Alpha_1
plot_title = ggplot2::ggtitle("Posterior distribution alpha_1", "with medians and 80% intervals")
mcmc_areas(results$draws("alpha_1"), prob = 0.8) + plot_title

plot_title = ggplot2::ggtitle("Traces for alpha_1")
mcmc_trace(results$draws("alpha_1")) + plot_title

# Lambda term
plot_title = ggplot2::ggtitle("Posterior distribution lambda_0", "with medians and 80% intervals")
mcmc_areas(results$draws("lambda_0"), prob = 0.8) + plot_title

plot_title = ggplot2::ggtitle("Traces for lambda_0")
mcmc_trace(results$draws("lambda_0")) + plot_title

plot_title = ggplot2::ggtitle("Posterior distribution lambda_1", "with medians and 80% intervals")
mcmc_areas(results$draws("lambda_1"), prob = 0.8) + plot_title

plot_title = ggplot2::ggtitle("Traces for lambda_1")
mcmc_trace(results$draws("lambda_1")) + plot_title

# plot_title = ggplot2::ggtitle("Posterior distribution lambda", "with medians and 80% intervals")
# mcmc_areas(results$draws("lambda"), prob = 0.8) + plot_title

# plot_title = ggplot2::ggtitle("Traces for lambda")
# mcmc_trace(results$draws("lambda")) + plot_title

#### Check-up correlations
cor(results$draws(variables = c("alpha_0"), format = "draws_df")$alpha_0,
	results$draws(variables = c("alpha_1"), format = "draws_df")$alpha_1)

cor(results$draws(variables = c("alpha_0"), format = "draws_df")$alpha_0,
	results$draws(variables = c("lambda"), format = "draws_df")$lambda)

cor(results$draws(variables = c("alpha_1"), format = "draws_df")$alpha_1,
	results$draws(variables = c("lambda"), format = "draws_df")$lambda)

#### Check-up plots with pointwise values (mean). _ptw for pointwise
dbh = trees[speciesName_fr == "Tilleul à grandes feuilles", dbh_in_mm]

# results$print(variables = c("alpha", "beta", "epsilon"), digits = 10)

alpha0_ptw = mean(results$draws(variables = c("alpha_0"), format = "draws_df")$alpha_0)
alpha1_ptw = mean(results$draws(variables = c("alpha_1"), format = "draws_df")$alpha_1)
# lambda_ptw = median(results$draws(variables = c("lambda"), format = "draws_df")$lambda)
lambda1_ptw = median(results$draws(variables = c("lambda_1"), format = "draws_df")$lambda_1)

fct_alpha = function(x)
	return (alpha0_ptw*x^alpha1_ptw)
	# return (alpha0_ptw*exp(-exp(-alpha1_ptw*x)))

fct_lambda = function(x)
	return (lambda1_ptw*x)
	# return (exp(-x/lambda1_ptw))

avg = function(x)
	return (fct_alpha(x)/(fct_alpha(x) - 1));
	# return (fct_lambda(x)/(fct_alpha(x) - 1) + 1);

skewness = function(x)
{
	l = fct_lambda(x)
	a = fct_alpha(x)
	return (2*l^3*a*(a + 1)/((a - 3)*(a - 2)*(a - 1)^3));
}

kurtosis = function(x)
{
	l = fct_lambda(x)
	a = fct_alpha(x)
	return (3*l^4*a*(3*a^2 + a + 2)/((a - 4)*(a - 3)*(a - 2)*(a - 1)^4));
}

n_rep = 100

dt = data.table(rep_dbh = rep(dbh, each = n_rep), sampled = numeric(n_rep * length(dbh)))

for (i in 1:length(dbh))
{
	alpha_i = fct_alpha(dbh[i]);
	# lambda_i = fct_lambda(dbh[i]);
	# dt[((i - 1)*n_rep + 1):(i*n_rep), sampled := rpareto2(n_rep, min = 1, shape = lambda_i, scale = alpha_i)]
	dt[((i - 1)*n_rep + 1):(i*n_rep), sampled := rpareto1(n_rep, shape = alpha_i, min = 1)]
}

plot(dbh, trees[speciesName_fr == "Tilleul à grandes feuilles", ratio], ylim = c(0.75, 1.5),
	pch = 19, xlab = "Diameter t0", ylab = "ratio diameter t1/diameter t0", col = "#323232AA")
curve(avg(x), col = "#CD212A", from = 0, to = 1500, lwd = 2, add = TRUE)
points(dt[, rep_dbh], dt[, sampled], pch = 19, col = "#4DA6FF02")

diam_class = seq(min(dbh) - 1, max(dbh) + 1, by = 10)
var_class = numeric(length(diam_class))
var_class_sim = numeric(length(diam_class))
for (i in 1:(length(diam_class) - 1))
{
	var_class[i] = var(trees[dbh_in_mm >= diam_class[i] & dbh_in_mm < diam_class[i + 1], ratio])
	var_class_sim[i] = var(dt[rep_dbh >= diam_class[i] & rep_dbh < diam_class[i + 1], sampled])
}

plot(diam_class, var_class, pch = 19, col = "#CD212A55")
points(diam_class, var_class_sim, pch = 19, col = "#4DA6FF55")


range(dt[, sampled])
hist(dt[, sampled])
hist(dt[sampled < 100, sampled])
100*dt[sampled < 100, .N]/dt[, .N]

plot(dt[, rep_dbh], dt[, sampled], pch = 19, col = "#4DA6FF")
curve(avg(x), col = "#CD212A", from = 0, to = 1500, lwd = 2)

index_min = which.min(dbh)
var_i = epsilon_ptw/dbh[index_min]^gamma_ptw;
avg_i = fct(dbh[index_min]);
scale_i = (1 - avg_i)*(1 - 2*avg_i + avg_i^2 + var_i)/(1 - 2*avg_i + avg_i^2 - var_i); # lambda in stan
shape_i = 2*var_i/(-1 + 2*avg_i + avg_i^2 + var_i); # alpha in stan
sim = rpareto2(n_rep, min = 1, shape = shape_i, scale = scale_i)

plot(rep(dbh[index_min], n_rep), sim)
points(rep(dbh[index_min], 61), trees[speciesName_fr == "Tilleul à grandes feuilles" & dbh_in_mm == min(dbh), ratio], col = "red")

#?#################################################################
#*###########    MODEL USING EITZEL 2013'S FORMULA    #############
#?#################################################################
#### Quick stan model on the ratio
model = cmdstan_model("eitzel_dbh.stan")
maxIter = 4e3
n_chains = 3

# init_lambda = lapply(1:n_chains, init_fun, dbh = trees[speciesName_fr == "Tilleul à grandes feuilles", dbh_in_mm])

stanData = list(N = trees[speciesName_fr == "Tilleul à grandes feuilles", .N],
	dbh0 = trees[speciesName_fr == "Tilleul à grandes feuilles", dbh_in_mm],
	dbh1 = trees[speciesName_fr == "Tilleul à grandes feuilles", dbh_5_in_mm])

results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 200, chains = n_chains,
	iter_warmup = maxIter/2, iter_sampling = maxIter/2, save_warmup = TRUE, # init = init_lambda,
	max_treedepth = 11, adapt_delta = 0.95)

results$cmdstan_diagnose()

# Alpha
plot_title = ggplot2::ggtitle("Posterior distribution alpha", "with medians and 80% intervals")
mcmc_areas(results$draws("alpha"), prob = 0.8) + plot_title

plot_title = ggplot2::ggtitle("Traces for alpha")
mcmc_trace(results$draws("alpha")) + plot_title

# Beta
plot_title = ggplot2::ggtitle("Posterior distribution beta", "with medians and 80% intervals")
mcmc_areas(results$draws("beta"), prob = 0.8) + plot_title

plot_title = ggplot2::ggtitle("Traces for beta")
mcmc_trace(results$draws("beta")) + plot_title

#### Check-up plots with pointwise values (mean). _ptw for pointwise
dbh = trees[speciesName_fr == "Tilleul à grandes feuilles", dbh_in_mm]

alpha_ptw = mean(results$draws(variables = c("alpha"), format = "draws_df")$alpha)
beta_ptw = mean(results$draws(variables = c("beta"), format = "draws_df")$beta)
sigma_eps_ptw = mean(results$draws(variables = c("sigma_eps"), format = "draws_df")$sigma_eps)

avg = function(x)
	return (alpha_ptw + beta_ptw*x);

n_rep = 100

dt = data.table(rep_dbh = rep(dbh, each = n_rep), sampled = numeric(n_rep * length(dbh)))

for (i in 1:length(dbh))
{
	avg_i = avg(dbh[i]);
	dt[((i - 1)*n_rep + 1):(i*n_rep), sampled := rnorm(n = n_rep, mean = avg_i, sd = sigma_eps_ptw)]
}

dt[, ratio := sampled/rep_dbh]

plot(dbh, trees[speciesName_fr == "Tilleul à grandes feuilles", dbh_5_in_mm], # ylim = c(0.75, 1.5),
	pch = 19, xlab = "Diameter t0", ylab = "ratio diameter t1/diameter t0", col = "#323232AA")
curve(avg(x), col = "#CD212A", from = 0, to = 1500, lwd = 2, add = TRUE)
points(dt[, rep_dbh], dt[, sampled], pch = 19, col = "#4DA6FF01")


plot(dbh, trees[speciesName_fr == "Tilleul à grandes feuilles", ratio], ylim = c(0.75, 1.5),
	pch = 19, xlab = "Diameter t0", ylab = "ratio diameter t1/diameter t0", col = "#323232AA")
curve(avg(x), col = "#CD212A", from = 0, to = 1500, lwd = 2, add = TRUE)
points(dt[, rep_dbh], dt[, ratio], pch = 19, col = "#4DA6FF02")

diam_class = seq(min(dbh) - 1, max(dbh) + 1, by = 10)
var_class = numeric(length(diam_class))
var_class_sim = numeric(length(diam_class))
for (i in 1:(length(diam_class) - 1))
{
	var_class[i] = var(trees[dbh_in_mm >= diam_class[i] & dbh_in_mm < diam_class[i + 1], ratio])
	var_class_sim[i] = var(dt[rep_dbh >= diam_class[i] & rep_dbh < diam_class[i + 1], sampled])
}

plot(diam_class, var_class, pch = 19, col = "#CD212A55")
points(diam_class, var_class_sim, pch = 19, col = "#4DA6FF55")

#! ----------------------------------------------------- !#
#! ------------------ CRASH TEST ZONE ------------------ !#
#! ----------------------------------------------------- !#
#### Petite disgression pour la variance
## Hint on data
diam_class = seq(min(dbh) - 1, max(dbh) + 1, by = 10)
var_class = numeric(length(diam_class))
for (i in 1:(length(diam_class) - 1))
	var_class[i] = var(trees[dbh_in_mm >= diam_class[i] & dbh_in_mm < diam_class[i + 1], ratio])

plot(diam_class, var_class)

## Model
model = cmdstan_model("var_class.stan")
maxIter = 4e3
n_chains = 3

stanData = list(N = length(var_class),
	var_class = var_class,
	diam_class = diam_class)

results = model$sample(data = stanData, parallel_chains = n_chains, refresh = 200, chains = n_chains,
	iter_warmup = maxIter/2, iter_sampling = maxIter/2, save_warmup = TRUE,
	max_treedepth = 12, adapt_delta = 0.95)

results$cmdstan_diagnose()

## A term
plot_title = ggplot2::ggtitle("Posterior distribution a", "with medians and 80% intervals")
mcmc_areas(results$draws("a"), prob = 0.8) + plot_title

# B term
plot_title = ggplot2::ggtitle("Posterior distribution b", "with medians and 80% intervals")
mcmc_areas(results$draws("b"), prob = 0.8) + plot_title

## Sigma term
plot_title = ggplot2::ggtitle("Posterior distribution sigma", "with medians and 80% intervals")
mcmc_areas(results$draws("sigma"), prob = 0.8) + plot_title

## Plot
a_ptw = mean(results$draws(variables = c("a"), format = "draws_df")$a)
b_ptw = mean(results$draws(variables = c("b"), format = "draws_df")$b)
sigma_ptw = mean(results$draws(variables = c("sigma"), format = "draws_df")$sigma)

fct = function(x)
	return (a_ptw/x^b_ptw)
	# return (a_ptw * exp(-b_ptw * x))

plot(diam_class, var_class)
curve(fct, from = min(diam_class), to = max(diam_class), lwd = 2, col = "#CD212A", add = TRUE)

fct_avg = function(x)
	return (a_ptw/x^b_ptw);

n_rep = 100
dt = data.table(rep_dbh = rep(dbh, each = n_rep), sampled = numeric(n_rep * length(dbh)))

for (i in 1:length(dbh))
{
	avg_i = fct_avg(dbh[i]);
	dt[((i - 1)*n_rep + 1):(i*n_rep), sampled := rnorm(n_rep, mean = avg_i, sd = sigma_ptw)]
}

plot(diam_class, var_class)
curve(fct, from = min(diam_class), to = max(diam_class), lwd = 2, col = "#CD212A", add = TRUE)
points(dt[, rep_dbh], dt[, sampled], pch = 19, col = "#4DA6FF03")

#### Test
# A = 1.54
# B = 2.9
# C = 1
# ll = rpareto2(n = 1e6, min = C, shape = B, scale = A)

# mean(ll)
# m = A/(B - 1) + C

# var(ll)
# v = A^2*B/((B - 2)*(B - 1)^2)

# lambda = (1 - m)*(1 - 2*m + m^2 + v)/(1 - 2*m + m^2 - v) # A
# alpha = 2*v/(-1 + 2*m - m^2 + v) # B

fct1 = function(x, alpha, theta, mu)
	return (alpha/(theta*(1 + (x - mu)/theta)^(alpha + 1)));

fct2 = function(x, alpha, lambda, mu) # http://www.naturalspublishing.com/files/published/73t3w33ckp06vv.pdf
	return (alpha/lambda * (1 + (x - mu)/lambda)^(-alpha - 1));

curve(fct2(x, 1, 1, 3), col = "#CD212A", from = 3, to = 15, lwd = 2, ylab = "Pareto type II")
curve(fct2(x, 1, 2, 3), col = "#FFA500", from = 3, to = 15, lwd = 2, add = TRUE)
curve(fct2(x, 1, 5, 3), col = "#E8A798", from = 3, to = 15, lwd = 2, add = TRUE)
curve(fct2(x, 2, 1, 3), col = "#0072B5", from = 3, to = 15, lwd = 2, add = TRUE)
curve(fct2(x, 5, 1, 3), col = "#34568B", from = 3, to = 15, lwd = 2, add = TRUE)
curve(fct2(x, 5, 25, 3), col = "#34568B", from = 3, to = 15, lwd = 2, add = TRUE)

qq = runif(n = 20, min = 2, max = 50)

max(abs(fct1(qq, 1.25, 0.5, 1) - fct2(qq, 1.25, 0.5, 1)))



lambda = function(avg, var)
	return ((1 - avg)*(1 - 2*avg + avg^2 + var)/(1 - 2*avg + avg^2 - var));

alpha = function(avg, var)
	return (2*var/(-1 + 2*avg - avg^2 + var));

min_x = 1.1
max_x = 5
min_y = 0
max_y = 1

x = seq(min_x, max_x, length = 10)
y = seq(min_y, max_y, length = 20)

z_l = outer(x, y, lambda)
sum(is.nan(z_l))
z_l_axis = unique(sort(c(0, seq(min(z_l), max(z_l), length.out = 5))))

mat_lambda = persp(x, y, z_l,
	theta = 40, phi = 30,
	ticktype = "detailed", box = FALSE, axes = FALSE, 
    mar = c(10, 10, 10, 20), expand = 0.5,
	col = "springgreen", shade = 0.5)

lines(trans3d(x, min_y, min(z_l), mat_lambda) , col = "black")
lines(trans3d(max_x, y, min(z_l), mat_lambda) , col = "black")
lines(trans3d(min_x, min_y, z_l_axis, mat_lambda) , col = "black")

tick_start = trans3d(x, min_y, min(z_l_axis), mat_lambda)
tick_end = trans3d(x, (min_y - 0.20), min(z_l_axis), mat_lambda)
segments(tick_start$x, tick_start$y, tick_end$x, tick_end$y)

tick_start = trans3d(max_x, y, min(z_l_axis), mat_lambda)
tick_end = trans3d(max_x + 0.20, y, min(z_l_axis), mat_lambda)
segments(tick_start$x, tick_start$y, tick_end$x, tick_end$y)

tick_start = trans3d(min_x, min_y, z_l_axis, mat_lambda)
tick_end = trans3d(min_x, (min_y - 0.1), z_l_axis, mat_lambda)
segments(tick_start$x, tick_start$y, tick_end$x, tick_end$y)

labels = as.character(ceiling(z_l_axis))
label_pos = trans3d(min_x, (min_y - 0.5), z_l_axis, mat_lambda)
text(label_pos$x, label_pos$y, labels = labels, adj = c(1, NA), cex = 0.75)

z_a = outer(x, y, alpha)



z_a = outer(x, y, alpha)

persp(x, y, z_a,
	xlab = "avg", ylab = "var", zlab = "alpha",
	theta = 0, phi = -10,
	col = "springgreen", shade = 0.5)

#### Pareto type 1
mean_par1 = function(alpha)
	return (alpha/(alpha - 1))

var_par1 = function(alpha)
	return (alpha/((alpha - 2)*(alpha - 1)^2))

curve(mean_par1, from = 1, to = 50)
curve(var_par1, from = 2, to = 50)
