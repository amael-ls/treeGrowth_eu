
#### Aim of prog: Plot different visualisations of growth versus climate and competition

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)

#### Tool functions
## Growth function
growth_fct = function(dbh, precipitation, temperature, basalArea, params, scaled_dbh = FALSE, scaled_clim = FALSE, ...)
{
	providedArgs = list(...)
	providedArgs_names = names(providedArgs)
	scaling = 1

	if (!scaled_dbh)
	{
		if (!("sd_dbh" %in% providedArgs_names))
			stop("You need to provide sd_dbh in order to standardise dbh")
		
		sd_dbh = providedArgs[["sd_dbh"]]
		scaling = sd_dbh
		dbh = dbh/sd_dbh
	}

	if (!scaled_clim)
	{
		if (!(c("pr_mu", "pr_sd", "tas_mu", "tas_sd") %in% providedArgs_names))
			stop("You need to provide pr_mu, pr_sd, tas_mu, and tas_sd in order to standardise dbh")
		
		temperature = (temperature - tas_mu)/tas_sd
		precipitation = (precipitation - pr_mu)/pr_sd
	}

	potentialGrowth = params["potentialGrowth"]
	dbh_slope = params["dbh_slope"]

	pr_slope = params["pr_slope"]
	pr_slope2 = params["pr_slope2"]
	tas_slope = params["tas_slope"]
	tas_slope2 = params["tas_slope2"]

	competition_slope = params["competition_slope"]

	return(scaling*exp(potentialGrowth + dbh_slope*dbh + pr_slope*precipitation + pr_slope2*precipitation^2 +
			tas_slope*temperature + tas_slope2*temperature^2 + competition_slope*basalRea))
}
