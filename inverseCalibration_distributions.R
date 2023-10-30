
#### Aim of prog: Compute different moments from mean and variance for few distributions
## Comments:
# It is easy to prove that the skewness of the gamma is always less than the skewness of the lognormal.
# Skewness(gamma) = 2*sd/mean
# Skewness(lognormal) = (var/mean^2 + 3)*sd/mean
# And therefore the difference (gamma - lognormal) is - sd/mean*(var/mean^2 + 1) < 0

library(extraDistr)
library(rootSolve)
library(nakagami)
library(moments)

inverseCalibration = function(fun, ...)
{
	# Check if distribution is included
	if (!isTRUE(all.equal(fun, dchisq)) &&
		!isTRUE(all.equal(fun, dgamma)) &&
		!isTRUE(all.equal(fun, dlnorm)) &&
		!isTRUE(all.equal(fun, dnaka)) &&
		!isTRUE(all.equal(fun, dnorm)) &&
		!isTRUE(all.equal(fun, dwald)) &&
		!isTRUE(all.equal(fun, dweibull)))
		stop("This function only accepts dchisq, dgamma, dlnorm, dnaka, dnorm, dwald or dweibull as priors")
	
	# Get list of arguments
	providedArgs = list(...)
	ls_names = names(providedArgs)

	# Define output
	output = list(mean = NA, sd = NA, var = NA, skewness = NA, arg1 = NA, arg2 = NA)

	# Chi-square
	if (isTRUE(all.equal(fun, dchisq))) # Checked and validated computation
	{
		if ((!all(c("mean", "var") %in% names(providedArgs))) && (!all(c("df", "ncp") %in% names(providedArgs))))
			stop("You must provide either mean and var or df and ncp for dchisq")
		
		if (all(c("mean", "var") %in% names(providedArgs)))
		{
			temp1 = providedArgs[["mean"]]
			temp2 = providedArgs[["var"]]

			arg1 = 2*temp1 - temp2/2 # degrees of freedom
			arg2 = temp2/2 - temp1 # non-central parameter
		} else {
			arg1 = providedArgs[["df"]]
			arg2 = providedArgs[["ncp"]]
		}
		output[["mean"]] = arg1 + arg2
		output[["var"]] = 2*arg1 + 4*arg2
		output[["sd"]] = sqrt(output[["var"]])
		output[["skewness"]] = sqrt(8)*(arg1 + 3*arg2)/(arg1 + 2*arg2)^1.3
		output[["arg1"]] = arg1
		output[["arg2"]] = arg2
	}

	# Gamma
	if (isTRUE(all.equal(fun, dgamma))) # Checked and validated computation
	{
		if ((!all(c("mean", "var") %in% names(providedArgs))) && (!all(c("shape", "rate") %in% names(providedArgs))))
			stop("You must provide either mean and var or shape and rate for dgamma")
		
		if (all(c("mean", "var") %in% names(providedArgs)))
		{
			temp1 = providedArgs[["mean"]]
			temp2 = providedArgs[["var"]]

			arg1 = temp1^2/temp2 # shape
			arg2 = temp1/temp2 # rate
		} else {
			arg1 = providedArgs[["shape"]]
			arg2 = providedArgs[["rate"]]
		}
		output[["mean"]] = arg1/arg2
		output[["var"]] = arg1/arg2^2
		output[["sd"]] = sqrt(output[["var"]])
		output[["skewness"]] = 2/sqrt(arg1) # = 2*sd/mean
		output[["arg1"]] = arg1
		output[["arg2"]] = arg2
	}

	# Lognormal
	if (isTRUE(all.equal(fun, dlnorm))) # Checked and validated computation
	{
		if ((!all(c("mean", "sd") %in% names(providedArgs))) && (!all(c("meanlog", "sdlog") %in% names(providedArgs))))
			stop("You must provide mean and sd or meanlog and sdlog for dlnorm")
		
		if (all(c("mean", "sd") %in% names(providedArgs)))
		{
			dlnorm_mean = providedArgs[["mean"]]
			dlnorm_sd = providedArgs[["sd"]]

			arg1 = log(dlnorm_mean^2/sqrt(dlnorm_sd^2 + dlnorm_mean^2))
			arg2 = sqrt(log(dlnorm_sd^2/dlnorm_mean^2 + 1))
		} else {
			arg1 = providedArgs[["meanlog"]]
			arg2 = providedArgs[["sdlog"]]
		}

		m = exp(arg1)
		omega = exp(arg2^2) # For compactness

		output[["mean"]] = m*exp(arg2^2/2)
		output[["var"]] = m^2*omega*(omega - 1)
		output[["sd"]] = sqrt(output[["var"]])
		output[["skewness"]] = (omega + 2)*sqrt(omega - 1) # = (var/mean^2 + 3)*var/mean
		output[["arg1"]] = arg1
		output[["arg2"]] = arg2
	}

	# Nakagami
	if (isTRUE(all.equal(fun, nakagami::dnaka))) # Checked and validated computation (with a trick for arg2)
	{
		if ((!all(c("mean", "var") %in% names(providedArgs))) && (!all(c("shape", "scale") %in% names(providedArgs))))
			stop("You must provide either mean and var or shape and rate for dnaka")
		
		if (all(c("mean", "var") %in% names(providedArgs)))
		{
			temp1 = providedArgs[["mean"]]
			temp2 = providedArgs[["var"]]

			eqnarray = function(shape, ...)
			{
				providedArgs = list(...)
				meanArg = providedArgs[["meanArg"]]
				varArg = providedArgs[["varArg"]]

				F1 = meanArg^2*(shape*gamma(shape)^2/gamma(shape + 0.5)^2 - 1) - varArg

				return (c(F1 = F1))
			}

			if ("start" %in% ls_names)
			{
				start = providedArgs[["start"]]
			} else {
				start = 5
			}

			arg1 = multiroot(f = eqnarray, start = start, meanArg = temp1, varArg = temp2)$root # shape
			arg2 = temp1*sqrt(arg1)*gamma(arg1)/gamma(arg1 + 0.5) # scale
		} else {
			arg1 = providedArgs[["shape"]]
			arg2 = providedArgs[["scale"]]
		}
		output[["mean"]] = arg2*gamma(arg1 + 0.5)/(sqrt(arg1)*gamma(arg1))
		output[["var"]] = arg2^2*(1 - gamma(arg1 + 0.5)^2/(arg1 * gamma(arg1)^2))
		output[["sd"]] = sqrt(output[["var"]])
		output[["skewness"]] = arg2^3*(2*gamma(arg1 + 0.5)^3 + 0.5*(1 - 4*arg1)*gamma(arg1)^2*gamma(arg1 + 0.5))/(arg1^1.5*gamma(arg1)^3)
		output[["arg1"]] = arg1
		# Necessary to put a square here, as dnaka uses the Wikipedia paramtrisation, while I used another one for analytical calculus
		output[["arg2"]] = arg2^2
	}

	# Normal
	if (isTRUE(all.equal(fun, dnorm))) # Checked and validated computation
	{
		if (!all(c("mean", "sd") %in% names(providedArgs)))
			stop("You must provide mean and sd for dnorm")
		
		arg1 = providedArgs[["mean"]]
		arg2 = providedArgs[["sd"]]

		output[["mean"]] = arg1
		output[["sd"]] = arg2
		output[["var"]] = arg2^2
		output[["skewness"]] = 0
		output[["arg1"]] = arg1
		output[["arg2"]] = arg2
	}

	# # Pareto
	# if (isTRUE(all.equal(fun, dpareto))) # Checked and validated computation
	# {
	# 	if ((!all(c("mean", "var") %in% names(providedArgs))) & (!all(c("location", "scale") %in% names(providedArgs))))
	# 		stop("You must provide either mean and var or shape and scale for dwald")

	# 	if (all(c("mean", "var") %in% names(providedArgs)))
	# 	{
	# 		arg1 = providedArgs[["mean"]] # location
	# 		arg2 = arg1^3/providedArgs[["var"]] # scale
	# 	} else {
	# 		arg1 = providedArgs[["location"]]
	# 		arg2 = providedArgs[["scale"]]
	# 	}
		
	# 	output[["mean"]] = arg1
	# 	output[["var"]] = arg1^3/arg2
	# 	output[["sd"]] = sqrt(output[["var"]])
	# 	output[["skewness"]] = 3*sqrt(arg1/arg2)
	# 	output[["arg1"]] = arg1
	# 	output[["arg2"]] = arg2
	# }

	# Wald
	if (isTRUE(all.equal(fun, dwald))) # Checked and validated computation
	{
		if ((!all(c("mean", "var") %in% names(providedArgs))) && (!all(c("location", "scale") %in% names(providedArgs))))
			stop("You must provide either mean and var or shape and scale for dwald")

		if (all(c("mean", "var") %in% names(providedArgs)))
		{
			arg1 = providedArgs[["mean"]] # location
			arg2 = arg1^3/providedArgs[["var"]] # scale
		} else {
			arg1 = providedArgs[["location"]]
			arg2 = providedArgs[["scale"]]
		}
		
		output[["mean"]] = arg1
		output[["var"]] = arg1^3/arg2
		output[["sd"]] = sqrt(output[["var"]])
		output[["skewness"]] = 3*sqrt(arg1/arg2)
		output[["arg1"]] = arg1
		output[["arg2"]] = arg2
	}

	# Weibull
	if (isTRUE(all.equal(fun, dweibull))) # Checked and validated computation
	{
		if ((!all(c("mean", "var") %in% names(providedArgs))) && (!all(c("shape", "scale") %in% names(providedArgs))))
			stop("You must provide either mean and var or shape and scale for dweibull")

		# There is no analytical method. Check https://www.scionresearch.com/__data/assets/pdf_file/0003/59286/NZJFS1131981GARCIA304-306.pdf
		if (all(c("mean", "var") %in% names(providedArgs)))
		{
			temp1 = providedArgs[["mean"]]
			temp2 = providedArgs[["var"]]

			eqnarray = function(shape, ...)
			{
				providedArgs = list(...)
				meanArg = providedArgs[["meanArg"]]
				varArg = providedArgs[["varArg"]]

				F1 = gamma(1 + 2/shape)/gamma(1 + 1/shape)^2 - 1 - varArg/meanArg^2

				return (c(F1 = F1))
			}

			if ("start" %in% ls_names)
			{
				start = providedArgs[["start"]]
			} else {
				start = 5
			}

			arg1 = multiroot(f = eqnarray, start = start, meanArg = temp1, varArg = temp2)$root # shape
			arg2 = temp1/gamma(1 + 1/arg1) # scale
		} else {
			arg1 = providedArgs[["shape"]]
			arg2 = providedArgs[["scale"]]
		}
		
		output[["mean"]] = arg2 * gamma(1 + 1/arg1)
		output[["var"]] = arg2^2 * ( gamma(1 + 2/arg1) - gamma(1 + 1/arg1)^2 )
		output[["sd"]] = sqrt(output[["var"]])
		output[["skewness"]] = ( 2*gamma(1 + 1/arg1)^3 - 3*gamma(1 + 1/arg1)*gamma(1 + 2/arg1) + gamma(1 + 3/arg1) )/
			(gamma(1 + 2/arg1) - gamma(1 + 1/arg1)^2)^(3/2)
		output[["arg1"]] = arg1
		output[["arg2"]] = arg2
	}

	return (output)
}

aa = inverseCalibration(dchisq, mean = 4, var = 8)
qq = rchisq(n = 1e7, df = aa[["arg1"]], ncp = aa[["arg2"]])

mean(qq)
var(qq)
skewness(qq)
