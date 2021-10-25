
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
	int<lower = 1> n_precip; // Dimension of the climate vector
	int<lower = 1> n_obs; // Number of trees observations
	int<lower = 2> n_hiddenState; // Dimension of the state space vector
	int<lower = n_obs - 1, upper = n_obs - 1> n_children; // Number of children trees observations = n_obs - 1
	int<lower = 2> nbYearsPerIndiv; // Number of years for each individual

	// Indices
	int<lower = 1, upper = n_obs - 1> parents_index; // Index of each parent in the observed data
	int<lower = 1, upper = n_hiddenState - 1> parentsObs_index; // Corresponding index of observed parents in Y_generated
	int<lower = 2, upper = n_obs> children_index[n_children]; // Index of children in the observed data
	int<lower = 2, upper = n_hiddenState> childrenObs_index[n_children]; // Corresponding index of observed children in Y_generated
	int<lower = 1, upper = n_precip - 1> climate_index; // Index of the climate associated to the parent (should be 1)
	int<lower = 2, upper = n_hiddenState> not_parent_index[n_hiddenState - 1]; // Index of states without data

	// Initial state, used to start the states and compare with Y_generated[parentsObs_index] (i.e., a prior parameter)
	real Y_generated_0;

	// Observations
	vector[n_obs] Yobs;

	// Explanatory variable
	vector<lower = 0>[n_precip] precip; // Precipitations
}

transformed data {
	vector[n_precip] normalised_precip = precip; // (precip - mean(precip))/sd(precip); // Normalised precipitations
}

parameters {
	// Parameters
	real slopes_dbh; // Species-specific value
	real slopes_precip; // Species-specific value
	real quad_slopes_precip; // Species-specific value

	real<lower = 0.000001> processError; // Constrained by default
	real<lower = 0.000001> measureError; // Constrained by default

	vector[n_hiddenState] Y_generated; // State vector Y (hidden markov process), positive
}

// THE COMPLEX INDEXING HAS BEEN CHECKED (23rd August 2021)
model {
	// Declare variables
	int count = 0;
	real expected_true_dbh[n_hiddenState - 1];
	int k = 0;

	// Priors
	// target += normal_lpdf(intercepts | 2, 10);
	target += normal_lpdf(slopes_dbh | 1, 1); // gamma_lpdf(slopes_dbh | 1.0^2/100, 1.0/100); // Gives a mean = 1, and variance = 100, /!\ integer division rounds to integer!
	target += normal_lpdf(slopes_precip | 0, 10);
	target += normal_lpdf(quad_slopes_precip | 0, 10);

	// print("target 3:", target());

	target += gamma_lpdf(processError | 10.0^2/1000, 10.0/1000); // Gives a mean  of 10 and variance of 1000
	target += gamma_lpdf(measureError | 2.0^2/10, 2.0/10); // Gives a mean  of 10 and variance of 1000
	// print("target 4:", target());

	// Model // THE COMPLEX INDEXING HAS BEEN CHECKED (23rd August 2021)
	for (i in 1:1)
	{
		for (j in 2:nbYearsPerIndiv) // Loop for all years but the first (which is the parent of indiv i)
		{
			k = k + 1;
			expected_true_dbh[k] = slopes_dbh*Y_generated[count + j - 1] + slopes_precip*normalised_precip[climate_index + j - 2] +
				quad_slopes_precip*normalised_precip[climate_index + j - 2]^2;
		}
		count += nbYearsPerIndiv;
	}
	// Prior on initial hidden state: This is a diffuse initialisation
	target += normal_lpdf(Y_generated[parentsObs_index] | Y_generated_0, processError); // should be processError instead of 10

	// Process model
	target += normal_lpdf(Y_generated[not_parent_index] | expected_true_dbh, processError); // should be processError instead of 10
	
	// --- Observation model
	// Compare true (hidden/latent) parents with observed parents
	target += normal_lpdf(Yobs[parents_index] | Y_generated[parentsObs_index], measureError); // should be measureError instead of 4

	// Compare true (hidden/latent) children with observed children
	target += normal_lpdf(Yobs[children_index] | Y_generated[childrenObs_index], measureError); // should be measureError instead of 4
}
