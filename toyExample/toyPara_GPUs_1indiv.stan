
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
	int<lower = 1> n_obs; // Number of trees observations
	int<lower = 2> n_hiddenState; // Dimension of the state space vector
	int<lower = n_obs - 1, upper = n_obs - 1> n_children; // Number of children trees observations = n_obs - 1
	int<lower = 2> nbYearsPerIndiv; // Number of years for each individual

	// Indices
	int<lower = 1, upper = n_obs - 1> parents_index; // Index of each parent in the observed data
	int<lower = 1, upper = n_hiddenState - 1> parentsObs_index; // Corresponding index of observed parents in Y_generated
	int<lower = 2, upper = n_obs> children_index[n_children]; // Index of children in the observed data
	int<lower = 2, upper = n_hiddenState> childrenObs_index[n_children]; // Corresponding index of observed children in Y_generated
	int<lower = 2, upper = n_hiddenState> not_parent_index[n_hiddenState - 1]; // Index of states without data

	// Initial state, used to start the states and compare with Y_generated[parentsObs_index] (i.e., a prior parameter)
	real Y_generated_0;

	// Observations
	vector[n_obs] Yobs;
}

parameters {
	// Parameters
	real<lower = 0.000001> processError; // Constrained by default
	real<lower = 0.000001> measureError; // Constrained by default
	real slope_dbh; // Constrained by default

	vector[n_hiddenState] Y_generated; // State vector Y (hidden markov process), positive
}

// THE COMPLEX INDEXING HAS BEEN CHECKED (23rd August 2021)
model {
	// Priors
	target += gamma_lpdf(processError | 10.0^2/1000, 10.0/1000); // Gives a mean  of 10 and variance of 1000
	target += gamma_lpdf(measureError | 10.0^2/1000, 10.0/1000); // Gives a mean  of 10 and variance of 1000
	target += normal_lpdf(slope_dbh | 1, 1);

	// Model
	target += normal_lpdf(Y_generated[1] | Y_generated_0, processError);
	for (i in 2:nbYearsPerIndiv)
	{
		target += normal_lpdf(Y_generated[i] | slope_dbh*Y_generated[i - 1], processError);
	}
	// --- Observation model
	// Compare true (hidden/latent) parents with observed parents
	target += normal_lpdf(Yobs[parents_index] | Y_generated[parentsObs_index], measureError);

	// Compare true (hidden/latent) children with observed children
	target += normal_lpdf(Yobs[children_index] | Y_generated[childrenObs_index], measureError);
}
