
/*
	Comments:
		- Yobs is an array of reals, containing all the observations
		- I call "parent" the first measured dbh of each individual. The parents have to be treated separately
	
	Reminder:
		- The gamma distribution of stan uses, in this order, shape (alpha) and rate (beta). It can however be reparametrised
			using mean and var: alpha = mean^2/var, beta = mean/var
*/

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
	int<lower = 2, upper = n_obs> children_index[n_children]; // Index of children in the observed data
	int<lower = 2, upper = n_hiddenState> childrenObs_index[n_children]; // Corresponding index of observed children in Y_generated
	int<lower = 1, upper = n_precip - 1> climate_index[n_indiv]; // Index of the climate associated to each parent

	// Observations
	vector<lower = 0>[n_obs] Yobs;

	// Explanatory variable
	vector<lower = 0>[n_precip] precip; // Precipitations
}

transformed data {
	vector[n_precip] normalised_precip = (precip - mean(precip))/sd(precip); // Normalised precipitations
}

parameters {
	// Parameters
	real intercepts; // Species-specific value
	real slopes_dbh; // Species-specific value, Markov process coefficient
	real slopes_precip; // Species-specific value
	real quad_slopes_precip; // Species-specific value

	real<lower = 0.000001> processError; // Constrained by default
	real<lower = 0.000001> measureError; // Constrained by default

	vector<lower = 0>[n_hiddenState] Y_generated; // State vector Y (hidden markov process), positive
}

// transformed parameters {
// 	real<lower = 0> quad_slopes_precip_m[n_indiv];
// 	for (i in 1:n_indiv)
// 		quad_slopes_precip_m[i] = -quad_slopes_precip[i];
// }

model {
	// Declare variables
	int count = 0;
	real mean_gamma_ij;

	// Priors
	target += normal_lpdf(intercepts | 0, 1000);
	target += gamma_lpdf(slopes_dbh | 1.0^2/10, 1.0/10); // Gives a mean of 1, and variance of 10, /!\ integer division rounds to integer!
	target += normal_lpdf(slopes_precip | 0, 1000);
	target += normal_lpdf(quad_slopes_precip | 0, 1000);

	// print("target 3:", target());

	target += gamma_lpdf(processError | 0.1, 0.01); // Gives a mean  of 0.1/0.01 = 10 and variance of 0.1/0.01^2 = 1000
	// print("target 4:", target());

	// Model
	for (i in 1:n_indiv)
	{
		// --- Parent data (uncertainty on the initial measure, should be country-specific!)
		// Y_generated[count + 1] ~ normal(Yobs[parents_index[i]], 1); // Prior
		target += normal_lpdf(Yobs[parents_index[i]] | Y_generated[count + 1], measureError); // Can be moved after and done in para!!!

		for (j in 2:nbYearsPerIndiv[i]) // Loop for all years but the first (which is the parent of indiv i)
		{
			mean_gamma_ij = intercepts + slopes_dbh*Y_generated[count + j - 1] +
				slopes_precip*normalised_precip[climate_index[i] + j - 2] + // -1 (shift index) and -1 (previous year)
				quad_slopes_precip*normalised_precip[climate_index[i] + j - 2]^2; // -1 (shift index) and -1 (previous year)
			
			// if (mean_gamma_ij <= 0)
			// {
			// 	print("----------");
			// 	print("**! Indiv = ", i);
			// 	print("Intercept = ", intercepts[i]);
			// 	print("slope dbh = ", slopes_dbh[i]);
			// 	print("slope p = ", slopes_precip[i]);
			// 	print("slope p^2 = ", quad_slopes_precip[i]);
			// 	print("Ygen[", count + j - 1, "] = ", Y_generated[count + j - 1]);
			// 	print("norm_precip[", climate_index[i] + j - 2, "] = ", normalised_precip[climate_index[i] + j - 2]);
			// }
			// Y_generated[count + j] ~ gamma(mean_gamma_ij^2/processError, mean_gamma_ij/processError);
			Y_generated[count + j] ~ normal(mean_gamma_ij, processError); // Can be moved after and done in para!!!
		}
		count = count + nbYearsPerIndiv[i];
	}
	// for (i in 1:n_obs)
	// 	print("y[", i, "] = ", Y_generated[i]);

	// Compare generated children with obs children
	target += normal_lpdf(Yobs[children_index] | Y_generated[childrenObs_index], measureError); // Can be done in para!!!
}

// generated quantities {
// 	real predict[N];
// 	int count = 0;

// 	for (i in 1:n_indiv)
// 	{
// 		predict[count + 1] = Yobs[parents_index[i]]; // Initiate parent
// 		for (j in 2:nbYearsPerIndiv[i]) // Loop for all years but the first (which is the parent of indiv i)
// 			predict[count + j] = intercepts[i] + slopes_dbh[i]*predict[count + j - 1] + slopes_precip*precip[count + j - 1];
// 		count = count + nbYearsPerIndiv[i];
// 	}
// }
