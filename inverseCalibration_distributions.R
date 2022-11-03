
library(extraDistr)
library(rootSolve)
library(moments)

inverseCalibration = function(fun, ...)
{
	# Check if distribution is included
	if (!isTRUE(all.equal(fun, dgamma)) &
		!isTRUE(all.equal(fun, dlnorm)) &
		!isTRUE(all.equal(fun, dnorm)) &
		!isTRUE(all.equal(fun, dwald)) &
		!isTRUE(all.equal(fun, dweibull)))
		stop("This function only accepts dgamma, dlnorm, dnorm, dwald or dweibull as priors")
	
	# Get list of arguments
	providedArgs = list(...)
	ls_names = names(providedArgs)

	# Define output
	output = list(mean = NA, sd = NA, var = NA, skewness = NA, arg1 = NA, arg2 = NA)

	# Get parameters
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

	if (isTRUE(all.equal(fun, dlnorm))) # Checked and validated computation
	{
		if ((!all(c("mean", "sd") %in% names(providedArgs))) &  (!all(c("meanlog", "sdlog") %in% names(providedArgs))))
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
		output[["skewness"]] = (omega + 2)*sqrt(omega - 1)
		output[["arg1"]] = arg1
		output[["arg2"]] = arg2
	}

	if (isTRUE(all.equal(fun, dgamma))) # Checked and validated computation
	{
		if ((!all(c("mean", "var") %in% names(providedArgs))) & (!all(c("shape", "rate") %in% names(providedArgs))))
			stop("You must provide either mean and var or shape and rate for dgamma")
		
		if (all(c("mean", "var") %in% names(providedArgs)))
		{
			temp1 = providedArgs[["mean"]]
			temp2 = providedArgs[["var"]] # Squared because in this case it is a variance, not a std. dev.

			arg1 = temp1^2/temp2 # shape
			arg2 = temp1/temp2 # rate
		} else {
			arg1 = providedArgs[["shape"]]
			arg2 = providedArgs[["rate"]]
		}
		output[["mean"]] = arg1/arg2
		output[["var"]] = arg1/arg2^2
		output[["sd"]] = sqrt(output[["var"]])
		output[["skewness"]] = 2/sqrt(arg1)
		output[["arg1"]] = arg1
		output[["arg2"]] = arg2
	}

	if (isTRUE(all.equal(fun, dweibull))) # Checked and validated computation
	{
		if ((!all(c("mean", "var") %in% names(providedArgs))) & (!all(c("shape", "scale") %in% names(providedArgs))))
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

			arg1 = multiroot(f = eqnarray, start = 5, meanArg = temp1, varArg = temp2)$root # shape
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

	if (isTRUE(all.equal(fun, dwald))) # Checked and validated computation
	{
		if ((!all(c("mean", "var") %in% names(providedArgs))) & (!all(c("location", "scale") %in% names(providedArgs))))
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

	return (output)
}

