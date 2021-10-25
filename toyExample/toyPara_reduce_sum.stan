
/*
	Comments:
		- Yobs is an array of reals, containing all the observations
		- I call "parent" the first measured dbh of each individual. The parents have to be treated separately
	
	Reminder:
		- The gamma distribution of stan uses, in this order, shape (alpha) and rate (beta). It can however be reparametrised
			using mean and var: alpha = mean^2/var, beta = mean/var
*/


// Signature  : real reduce_sum(F f, T[] x, int grainsize, T1 s1, T2 s2, ...)
// Signature f: real f(T[] x_slice, int start, int end, T1 s1, T2 s2, ...) where T1, T2, ... are whatever parameters

functions {
	real partial_sum_lpdf(real[] slice_obs, int start, int end, real[] subset_Y_gen, real measureError)
		{
			return normal_lpdf(slice_obs | subset_Y_gen[start:end], measureError);
		}
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
	int<lower = 1, upper = n_hiddenState - 1> parentsObs_index[n_indiv]; // Corresponding index of observed parents in Y_generated
	int<lower = 2, upper = n_obs> children_index[n_children]; // Index of children in the observed data
	int<lower = 2, upper = n_hiddenState> childrenObs_index[n_children]; // Corresponding index of observed children in Y_generated
	int<lower = 1, upper = n_precip - 1> climate_index[n_indiv]; // Index of the climate associated to each parent
	int<lower = 2, upper = n_hiddenState> not_parent_index[n_hiddenState - n_indiv]; // Index of states without data

	// Observations
	real<lower = 0> Yobs[n_obs];

	// Explanatory variable
	vector<lower = 0>[n_precip] precip; // Precipitations

	// Parameter for parallel calculus
	int<lower = 1, upper = n_hiddenState> grainsize; // cf https://mc-stan.org/docs/2_27/stan-users-guide/reduce-sum.html
}

transformed data {
	vector[n_precip] normalised_precip = (precip - mean(precip))/sd(precip); // Normalised precipitations
}

parameters {
	// Parameters
	real intercepts; // Species-specific value
	real<lower = 0> slopes_dbh; // Species-specific value, Markov process coefficient
	real slopes_precip; // Species-specific value
	real quad_slopes_precip; // Species-specific value

	real<lower = 0.000001> processError; // Constrained by default
	real<lower = 0.000001> measureError; // Constrained by default

	real<lower = 0> Y_generated[n_hiddenState]; // State vector Y (hidden markov process), positive
}

model {
	// Declare variables
	int count = 0;
	real expected_true_dbh[n_hiddenState - n_indiv];
	int k = 0;

	// Priors
	target += normal_lpdf(intercepts | 0, 1000);
	target += gamma_lpdf(slopes_dbh | 1.0^2/10, 1.0/10); // Gives a mean of 1, and variance of 10, /!\ integer division rounds to integer!
	target += normal_lpdf(slopes_precip | 0, 1000);
	target += normal_lpdf(quad_slopes_precip | 0, 1000);

	// print("target 3:", target());

	target += gamma_lpdf(processError | 0.1, 0.01); // Gives a mean  of 0.1/0.01 = 10 and variance of 0.1/0.01^2 = 1000
	target += gamma_lpdf(measureError | 0.1, 0.01); // Gives a mean  of 0.1/0.01 = 10 and variance of 0.1/0.01^2 = 1000
	// print("target 4:", target());

	// Model
	for (i in 1:n_indiv)
	{
		for (j in 2:nbYearsPerIndiv[i]) // Loop for all years but the first (which is the parent of indiv i)
		{
			k = k + 1;
			expected_true_dbh[k] = intercepts + slopes_dbh*Y_generated[count + j - 1] +
				slopes_precip*normalised_precip[climate_index[i] + j - 2] + // -1 (shift index) and -1 (previous year)
				quad_slopes_precip*normalised_precip[climate_index[i] + j - 2]^2; // -1 (shift index) and -1 (previous year)
		}
		count = count + nbYearsPerIndiv[i];
	}
	// Prior on initial hidden state
	target += reduce_sum(partial_sum_lpdf, Y_generated[parentsObs_index], grainsize, Yobs[parents_index], 1);

	// Process model
	target += reduce_sum(partial_sum_lpdf, Y_generated[not_parent_index], grainsize, expected_true_dbh, processError);
	
	// --- Observation model
	// Compare true (hidden/latent) parents with observed parents
	target += reduce_sum(partial_sum_lpdf, Yobs[parents_index], grainsize, Y_generated[parentsObs_index], measureError);

	// Compare true (hidden/latent) children with observed children
	target += reduce_sum(partial_sum_lpdf, Yobs[children_index], grainsize, Y_generated[childrenObs_index], measureError);
}
