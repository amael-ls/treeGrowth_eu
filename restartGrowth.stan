
// https://discourse.mc-stan.org/t/looking-for-advice-on-fitting-a-time-series-model-every-single-day/19065


/*
	Comments:
		- Yobs is an vector, containing all dbh observations
		- I call "parent" the first measured dbh of each individual. The parents have to be treated separately
		- The observation error for the French data is different from the other countries because there have been corrections in the data.
			Therefore, there is no extrem error due to, for instance, typos.
		- Note that the dbh is standardised but not centred. The reason was that in a previous version of the model, it would have
			implied a division by 0. Since then, I do not have this problem anymore, but I kept the dbh non-centred.
	
	Reminder:
		- The gamma distribution of stan uses, in this order, shape (alpha) and rate (beta). It can however be reparametrised
			using mean and var: alpha = mean^2/var, beta = mean/var
*/

functions {
	// Function for growth. This returns the expected growth in mm, for 1 year.
	real growth(real dbh0, real precip, real temperature, real totalTreeWeight,
		real potentialGrowth, real dbh_slope, real pr_slope, real pr_slope2 , real tas_slope, real tas_slope2, real competition_slope)
	{
		/*
			It takes the following variables:
				- dbh0: The diameter
				- precip: The precipitations
				- temperature: The temperature
				- totalTreeWeight: The sum of all the tree weights
			
			It takes the following parameters:
				- potentialGrowth: The basic growth, when all the explanatory variables are set to zero
				- dbh_slope: The slope for dbh
				- pr_slope: The slope for precipitation
				- pr_slope2: The slope for precipitation (squared term)
				- tas_slope: The slope for temperature (tas)
				- tas_slope2: The slope for temperature (tas, squared term)
				- competition_slope: The slope for competition
		*/

		return (exp(potentialGrowth + dbh_slope*dbh0 + pr_slope*precip + pr_slope2*precip^2 +
			tas_slope*temperature + tas_slope2*temperature^2 + competition_slope*totalTreeWeight));
	}
}

data {
	// Number of data
	int<lower = 1> n_indiv; // Number of individuals
	int<lower = 1> n_climate; // Dimension of the climate vector
	int<lower = 1> n_obs; // Number of trees observations
	int<lower = 1> n_latentGrowth; // Dimension of the state space vector for latent growth
	int<lower = n_obs - n_indiv, upper = n_obs - n_indiv> n_children; // Number of children trees observations = n_obs - n_indiv
	array [n_indiv] int<lower = 2, upper = n_obs> nbYearsGrowth; // Number of years of growth for each individual

	// Indices
	array [n_indiv] int<lower = 1, upper = n_obs - 1> parents_index; // Index of each parent in the observed data
	array [n_children] int<lower = 2, upper = n_obs> children_index; // Index of children in the observed data
	array [n_indiv] int<lower = 1, upper = n_climate - 1> climate_index; // Index of the climate associated to each parent

	// Observations
	vector<lower = 0>[n_obs] Yobs;

	// sd_dbh is for the whole species and should not be more than 5% different from the sd of the subsample, namely sd(Yobs)
	real<lower = 0.95*sd(Yobs), upper = 1.05*sd(Yobs)> sd_dbh;

	// Explanatory variables
	vector<lower = 0>[n_climate] precip; // Precipitations
	real pr_mu; // To standardise the precipitations
	real<lower = 0> pr_sd; // To standardise the precipitations

	vector[n_climate] tas; // Temperature
	real tas_mu; // To standardise the temperature
	real<lower = 0> tas_sd; // To standardise the temperature

	/*
		The tree weight correspond to the weight of the tree per hectar. It accounts for:
			1. the area of the plot (disc with a radius of 6, 9 or 15 metres depending on the size class)
			2. the proximity to the boundary (rectification of the probability of drawing trees at the edge)
	*/
	vector<lower = 0>[n_indiv] totalTreeWeight; // Sum of the tree weights for a given plot at a given tieme

	// Posterior from first run
	real potentialGrowthPrior_mean;
	real<lower = 0> potentialGrowthPrior_sd;

	real dbh_slopePrior_mean;
	real<lower = 0> dbh_slopePrior_sd;

	real pr_slopePrior_mean;
	real<lower = 0> pr_slopePrior_sd;
	real pr_slope2Prior_mean;
	real<lower = 0> pr_slope2Prior_sd;

	real tas_slopePrior_mean;
	real<lower = 0> tas_slopePrior_sd;
	real tas_slope2Prior_mean;
	real<lower = 0> tas_slope2Prior_sd;

	real competition_slopePrior_mean;
	real<lower = 0> competition_slopePrior_sd;

	real processErrorPrior_mean;
	real<lower = 0> processErrorPrior_sd;
	
	real measureErrorPrior_mean;
	real<lower = 0> measureErrorPrior_sd;
}

transformed data {
	vector[n_obs] normalised_Yobs = Yobs/sd_dbh; // Normalised but NOT centred dbh
	vector[n_climate] normalised_precip = (precip - pr_mu)/pr_sd; // Normalised and centred precipitations
	vector[n_climate] normalised_tas = (tas - tas_mu)/tas_sd; // Normalised and centred temperatures
	vector[n_indiv] normalised_totalTreeWeight = (totalTreeWeight - mean(totalTreeWeight))/sd(totalTreeWeight);
	// Normalised and centred totalTreeWeight
}

parameters {
	// Parameters
	real potentialGrowth; // Growth when all the explanatory variables are set to 0
	real dbh_slope;

	real pr_slope;
	real pr_slope2;

	real tas_slope;
	real tas_slope2;

	real competition_slope;

	real<lower = 0.5/sd_dbh^2> processError;
	real<lower = 0.1/sqrt(12)*25.4/sd_dbh> measureError; // Constrained by default, see appendix D Eitzel for the calculus

	vector<lower = 0.1/sd_dbh, upper = 2000/sd_dbh>[n_indiv] latent_dbh_parents; // Real (and unobserved) first measurement dbh (parents)
	vector<lower = 0>[n_latentGrowth] latent_growth; // Real (and unobserved) yearly growth
}

model {
	// Declare variables
	int growth_counter = 1;
	real expected_growth;
	vector [n_obs - n_indiv] latent_dbh_children;

	// Priors, the 1.1*..._sd is to increase the prior's variance of the first fit about 10%
	target += normal_lpdf(potentialGrowth | potentialGrowthPrior_mean, 1.1*potentialGrowthPrior_sd);
	target += normal_lpdf(dbh_slope | dbh_slopePrior_mean, 1.1*dbh_slopePrior_sd);
	
	target += normal_lpdf(pr_slope | pr_slopePrior_mean, 1.1*pr_slopePrior_sd);
	target += normal_lpdf(pr_slope2 | pr_slope2Prior_mean, 1.1*pr_slope2Prior_sd);

	target += normal_lpdf(tas_slope | tas_slopePrior_mean, 1.1*tas_slopePrior_sd);
	target += normal_lpdf(tas_slope2 | tas_slope2Prior_mean, 1.1*tas_slope2Prior_sd);

	target += normal_lpdf(competition_slope | competition_slopePrior_mean, 1.1*competition_slopePrior_sd);

	/*
		A note on the two folowwing priors:
		1. Because the values provided here are already scaled, there is no need to scale them again, although sd(Yobs) from this run is not
			exactly equal to sd(Yobs) from the first run. I secured this by adding a warning in the R code if there is more than 5%
			difference between the sd of the whole data set and the sd used for the first run. Note that in the stan code itself I require
			that sd_dbh is within +/- 5% of sd(Yobs)

		2. Note that ***ErrorPrior_sd is a standard deviation, but the gamma distribution needs the variance. Therefore, it is squared!

	*/
	target += gamma_lpdf(processError | processErrorPrior_mean^2/processErrorPrior_sd^2, processErrorPrior_mean/processErrorPrior_sd^2);
	target += gamma_lpdf(measureError | measureErrorPrior_mean^2/measureErrorPrior_sd^2, measureErrorPrior_mean/measureErrorPrior_sd^2);

	// Model
	for (i in 1:n_indiv) // Loop over all the individuals
	{
		latent_dbh_children[i] = latent_dbh_parents[i];
		for (j in 1:nbYearsGrowth[i]) // Loop for all years but the first (which is the parent of indiv i)
		{
			// Process model
			expected_growth = growth(latent_dbh_children[i], normalised_precip[climate_index[i] + j - 1],
				normalised_tas[climate_index[i] + j - 1], normalised_totalTreeWeight[i],
				potentialGrowth, dbh_slope, pr_slope, pr_slope2, tas_slope, tas_slope2, competition_slope);
			target += gamma_lpdf(latent_growth[growth_counter] | expected_growth^2/processError, expected_growth/processError);

			// Dbh at time t + 1, only the last (i.e., child) dbh is recorded
			latent_dbh_children[i] += latent_growth[growth_counter];
			growth_counter += 1;
		}
	}
	
	// Prior on initial hidden state: This is a diffuse initialisation
	target += uniform_lpdf(latent_dbh_parents | 0.1/sd_dbh, 2000/sd_dbh); // Remember that the dbh is in mm and standardised
	
	// --- Observation model
	// Compare true (i.e., hidden or latent) parents with observed parents
	target += normal_lpdf(normalised_Yobs[parents_index] | latent_dbh_parents, measureError);

	// Compare true (i.e., hidden or latent) children with observed children
	target += normal_lpdf(normalised_Yobs[children_index] | latent_dbh_children, measureError);
}
