
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

	// Explanatory variable
	vector<lower = 0>[n_precip] precip; // Precipitations
	real climate_mu;
	real<lower = 0> climate_sd;
}

transformed data {
	vector[n_obs] normalised_Yobs = Yobs/sd(Yobs); // Normalised but NOT centred dbh
	vector[n_precip] normalised_precip = (precip - climate_mu)/climate_sd; // Normalised precipitations
	
	// Constraints
	real lowest_prec = (0.0 - climate_mu)/climate_sd; // Lowest precipitation on the normalised scale
}

parameters {
	// Parameters
	real<lower = 0> potentialMaxGrowth;
	real power_dbh; // Coefficient showing how growth is a (hopefully!) decreasing function of dbh
	real<lower = lowest_prec> optimal_precip; // Optimal precipitation niche value
	real<lower = 0> width_precip_niche; // Width of the niche along the precipitation axis

	real<lower = 0.01, upper = 10> processError; // Constrained by default, realistically not too small
	// real<lower = 0.0055, upper = 0.011> measureError; // Constrained by default, lower bound = 0.1/sqrt(12)*25.4/sd(dbh). See appendix D Eitzel for the calculus
	// real<lower = 0.0055> measureError; // TEST FRENCH DATA

	vector[n_hiddenState] latent_dbh; // Real (and unobserved) dbh
}

model {
	// Declare variables
	int count = 0;
	real expected_latent_dbh[n_hiddenState - n_indiv];
	real expected_latent_yearlyGrowth;
	int k = 0;

	// Priors
	target += gamma_lpdf(potentialMaxGrowth | 1.0^2/10, 1.0/10);
	target += normal_lpdf(power_dbh | 0, 5);
	target += normal_lpdf(optimal_precip | 0, 20);
	target += gamma_lpdf(width_precip_niche | 1.0^2/10000, 1.0/10000); // Gives a mean  of 1 and variance of 10000

	target += gamma_lpdf(processError | 1.0^2/100, 1.0/100); // Gives a mean  of 1 and variance of 100
	// target += uniform_lpdf(measureError | 0.0055, 0.011); // The upper bound means that there is at max an error of 23.35 mm on the circumference
	// target += gamma_lpdf(measureError | 0.07^2/0.1, 0.07/0.1); // TEST FRENCH DATA
	// target += normal_lpdf(measureError | 3.0/135.137, 1.0/135.137); // TEST 2 FRENCH DATA: Correspond to a dbh measurement error of 3 mm, sd(dbh) = 135.137

	// Model
	for (i in 1:n_indiv)
	{
		for (j in 2:nbYearsPerIndiv[i]) // Loop for all years but the first (which is the parent of indiv i)
		{
			k = k + 1;
			expected_latent_yearlyGrowth = potentialMaxGrowth *
				latent_dbh[count + j - 1]^power_dbh *
				exp(-(normalised_precip[climate_index[i] + j - 2] - optimal_precip)^2/width_precip_niche^2);

			expected_latent_dbh[k] = latent_dbh[count + j - 1] + expected_latent_yearlyGrowth;
		}
		count += nbYearsPerIndiv[i];
	}
	// Prior on initial hidden state: This is a diffuse initialisation
	target += normal_lpdf(latent_dbh[parentsObs_index] | 0, 1e6); // Diffuse initialisation

	// Process model
	target += normal_lpdf(latent_dbh[not_parent_index] | expected_latent_dbh, processError);
	
	// --- Observation model
	// Compare true (hidden/latent) parents with observed parents
	// target += normal_lpdf(Yobs[parents_index] | latent_dbh[parentsObs_index], 0.795); // 2.5/pi
	target += normal_lpdf(normalised_Yobs[parents_index] | latent_dbh[parentsObs_index], 1.0/135.137); // measureError, 2.5/pi/sd(dbh) (WHY THIS?)

	// Compare true (hidden/latent) children with observed children
	// target += normal_lpdf(Yobs[children_index] | latent_dbh[childrenObs_index], 0.795); // 2.5/pi
	target += normal_lpdf(normalised_Yobs[children_index] | latent_dbh[childrenObs_index], 1.0/135.137); // measureError, 2.5/pi/sd(dbh) (WHY THIS?)
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
