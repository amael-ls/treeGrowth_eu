
/*
	Comments:
		- Yobs is an array of reals, containing all the observations
		- I call "parent" the first measured dbh of each individual. The parents have to be treated separately
		- The observation error for the French data is fixed to 2.5 mm on the circumference. This gives an error of 2.5/pi mm on the dbh
			Note that this error must then be expressed on the normalised scale: divide by sd(dbh)
		- This version of the growth model follows a discussion with Florian and Lisa
		- Note that I did not normalise dbh as I use a power function in the posterior distrib. This prevent a division by 0
	
	Reminder:
		- The gamma distribution of stan uses, in this order, shape (alpha) and rate (beta). It can however be reparametrised
			using mean and var: alpha = mean^2/var, beta = mean/var

	Information on data:
		mu(dbh) = 219.2785
		sd(dbh) = 135.137
*/

// PROCESS ERROR ON GROWTH RATHER THAN DBH

functions {
	// Function for growth. This returns the expected growth in mm, for 1 year.
	real growth(real dbh0, real precip, /* real temp, real total_TA, */
		real potentialGrowth, real dbh_slope, real pr_slope, real pr_slope2 /*, real temp_slope, real temp_slope2, real trunkArea_slope*/)
	{
		/*
			It takes the following variables:
				- dbh0: The diameter
				- precip: The precipitations
				- temp: The temperature
				- total_TA: The total trunk area
			
			It takes the following parameters:
				- potentialGrowth: The basic growth, when all the explanatory variables are set to zero
				- dbh_slope: The slope for dbh
				- pr_slope: The slope for precipitation
				- pr_slope2: The slope for precipitation (squared term)
				- temp_slope: The slope for precipitation
				- temp_slope2: The slope for precipitation (squared term)
				- trunkArea_slope: The slope for competition
		*/

		return (exp(potentialGrowth + dbh_slope*dbh0 + pr_slope*precip + pr_slope2*precip^2));
	}

	// real processError_fct(real dbh, real alpha, real beta)
	// {
	// 	return (exp(alpha + beta*dbh));
	// }
}

data {
	// Number of data
	int<lower = 1> n_indiv; // Number of individuals
	int<lower = 1> n_precip; // Dimension of the climate vector
	int<lower = 1> n_obs; // Number of trees observations
	int<lower = 2> n_hiddenState; // Dimension of the state space vector
	int<lower = n_obs - n_indiv, upper = n_obs - n_indiv> n_children; // Number of children trees observations = n_obs - n_indiv
	int<lower = 2, upper = n_obs> nbYearsPerIndiv[n_indiv]; // Number of years for each individual

	// Indices
	int<lower = 1, upper = n_obs - 1> parents_index[n_indiv]; // Index of each parent in the observed data
	int<lower = 1, upper = n_hiddenState - 1> parentsObs_index[n_indiv]; // Corresponding index of observed parents in latentState
	int<lower = 2, upper = n_obs> children_index[n_children]; // Index of children in the observed data
	int<lower = 2, upper = n_hiddenState> childrenObs_index[n_children]; // Corresponding index of observed children in latentState
	int<lower = 1, upper = n_precip - 1> climate_index[n_indiv]; // Index of the climate associated to each parent
	int<lower = 2, upper = n_hiddenState> not_parent_index[n_hiddenState - n_indiv]; // Index of states without data

	// Observations
	vector<lower = 0>[n_obs] Yobs;
	vector<lower = 0>[n_hiddenState] latent_dbh; // The answer is provided

	// Explanatory variables
	vector<lower = 0>[n_precip] precip; // Precipitations
	real pr_mu; // To standardise the precipitations
	real<lower = 0> pr_sd; // To standardise the precipitations

	vector<lower = 0>[n_obs] totalTrunkArea;
}

transformed data {
	vector[n_obs] normalised_Yobs = Yobs/sd(Yobs); // Normalised but NOT centred dbh
	vector[n_precip] normalised_precip = (precip - pr_mu)/pr_sd; // Normalised and centered precipitations
}

parameters {
	// Parameters
	real potentialGrowth; // Growth when all the explanatory variables are set to 0
	real dbh_slope;
	
	// real comp;

	real pr_slope;
	real pr_slope2;

	// real<lower = 0> processError; // Constrained by default, realistically not too small
	// real<lower = 0.1/sqrt(12)*25.4/sd(Yobs)> measureError; // Constrained by default, see appendix D Eitzel for the calculus

	// vector<lower = 0>[n_hiddenState] latent_dbh; // Real (and unobserved) dbh
}

model {
	// Declare variables
	real expected_latent_dbh[n_hiddenState - n_indiv];
	real measureError = 3.0/sd(Yobs); // measureError is a standard deviation
	real processError = 5.68/sd(Yobs)^2; // processError is a variance
	int count = 0;
	int k = 0;

	// Priors
	target += normal_lpdf(potentialGrowth | 0, 100);
	target += normal_lpdf(dbh_slope | 0, 5);
	
	// target += normal_lpdf(comp | 0, 5);

	target += normal_lpdf(pr_slope | 0, 5);
	target += normal_lpdf(pr_slope2 | 0, 5);

	// target += gamma_lpdf(processError | 1.0^2/10000, 1.0/10000); // Gives a mean  of 1 and variance of 10000
	// target += normal_lpdf(measureError | 3.0/sd(Yobs), 0.25/sd(Yobs)); // Correspond to a dbh measurement error of 3 mm, standardised

	// Model
	for (i in 1:n_indiv) // Loop over all the individuals
	{
		for (j in 2:nbYearsPerIndiv[i]) // Loop for all years but the first (which is the parent of indiv i)
		{
			k = k + 1;
			expected_latent_dbh[k] = latent_dbh[count + j - 1] + growth(latent_dbh[count + j - 1], normalised_precip[climate_index[i] + j - 2],
				potentialGrowth, dbh_slope, pr_slope, pr_slope2);
			// j - 2 comes from (j - 1) - 1. The first '-1' is to compensate that j starts at 2, the second '-1' is to match the latent_dbh
		}
		count += nbYearsPerIndiv[i];
	}
	// Prior on initial hidden state: This is a diffuse initialisation
	// target += gamma_lpdf(latent_dbh[parentsObs_index] | 1.582318^2/10.0, 1.582318/10.0); // Remember that the dbh is standardised

	// Process model
	// target += normal_lpdf(latent_dbh[not_parent_index] | expected_latent_dbh, processError);
	for (i in 1:(n_hiddenState - n_indiv))
		target += gamma_lpdf(latent_dbh[not_parent_index[i]] | expected_latent_dbh[i]^2/processError, expected_latent_dbh[i]/processError);
	
	// --- Observation model
	// Compare true (hidden/latent) parents with observed parents
	// target += normal_lpdf(normalised_Yobs[parents_index] | latent_dbh[parentsObs_index], measureError);

	// Compare true (hidden/latent) children with observed children
	// target += normal_lpdf(normalised_Yobs[children_index] | latent_dbh[childrenObs_index], measureError);
}

// Model with stuff written by Florian
// model {
// 	// Declare variables
// 	int count = 0;
// 	real expected_latent_dbh[n_hiddenState - n_indiv];
// 	int k = 0;

// 	// Priors
// 	target += normal_lpdf(intercepts | 0, 100);
// 	target += normal_lpdf(slopes_dbh | 0, 100);
// 	target += normal_lpdf(slopes_precip | 0, 100);
// 	target += normal_lpdf(quad_slopes_precip | 0, 100);

// 	// print("target 3:", target());

// 	target += gamma_lpdf(processError | 10.0^2/10000, 10.0/10000); // Gives a mean  of 10 and variance of 10000
// 	// target += gamma_lpdf(measureError | 1.0^2/1000, 1.0/1000); // Gives a mean  of 10 and variance of 1000

// 	// Model // THE COMPLEX INDEXING HAS BEEN CHECKED (23rd August 2021)
// 	for (i in 1:n_indiv)
// 	{
// 		for (j in 2:nbYearsPerIndiv[i]) // Loop for all years but the first (which is the parent of indiv i)
// 		{
// 			k = k + 1;
//			expectedGrowth[k] = exp(baseGrowth + growthL *latentState[count + j - 1] +
//				slopes_precip*normalised_precip[climate_index[i] + j - 2] + // -1 (shift index) and -1 (previous year)
//				quad_slopes_precip*normalised_precip[climate_index[i] + j - 2]^2); // -1 (shift index) and -1 (previous year)

// 			latentState[count + j] = latentState[count + j - 1] + realGrowth
// 		}
// 		count += nbYearsPerIndiv[i];
// 	}
// 	// Prior on initial hidden state: This is a diffuse initialisation
// 	target += normal_lpdf(latentState[parentsObs_index] | 0, 1e6); // Diffuse initialisation
// 	// target += normal_lpdf(latentState[parentsObs_index] | initialParents, 10); // Not too diffused initialisation

// 	// Process model
// 	target += gamma_lpdf(realGrowth | expectedGrowth, processError);
	
// 	// --- Observation model
// 	// Compare true (hidden/latent) parents with observed parents
// 	target += normal_lpdf(Yobs_normalised[parents_index] | latentState[parentsObs_index], 0.006); // (2.5/pi)/sd(dbh)

// 	// Compare true (hidden/latent) children with observed children
// 	target += normal_lpdf(Yobs_normalised[children_index] | latentState[childrenObs_index], 0.006); // (2.5/pi)/sd(dbh)
// }
