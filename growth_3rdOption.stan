
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
	int<lower = 2, upper = n_obs> children_index[n_children]; // Index of children in the observed data
	int<lower = 1, upper = n_precip - 1> climate_index[n_indiv]; // Index of the climate associated to each parent

	// Observations
	vector<lower = 0>[n_obs] Yobs;

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

	real pr_slope;
	real pr_slope2;

	// real<lower = 0> processError; // Constrained by default, realistically not too small
	real<lower = 0.1/sqrt(12)*25.4/sd(Yobs)> measureError; // Constrained by default, see appendix D Eitzel for the calculus

	vector<lower = 0>[n_indiv] latent_dbh_parents; // Real (and unobserved) parents dbh
	vector<lower = 0>[n_hiddenState - n_indiv] latent_growth; // Real (and unobserved) yearly growth
}

model {
	// Declare variables
	real expected_growth;
	real latent_dbh_children[n_indiv];
	real current_dbh;
	real processError = 5.68/sd(Yobs)^2; // processError is a variance
	int k = 0;

	// Priors
	target += normal_lpdf(potentialGrowth | 0, 100);
	target += normal_lpdf(dbh_slope | 0, 5);
	
	target += normal_lpdf(pr_slope | 0, 5);
	target += normal_lpdf(pr_slope2 | 0, 5);

	// target += gamma_lpdf(processError | 1.0^2/10000, 1.0/10000); // Gives a mean  of 1 and variance of 10000
	target += normal_lpdf(measureError | sqrt(1.0)/sd(Yobs), 0.15/sd(Yobs)); // Correspond to a dbh measurement error of 3 mm, standardised

	// Model
	for (i in 1:n_indiv) // Loop over all the individuals
	{
		current_dbh = latent_dbh_parents[i];
		latent_dbh_children[i] = latent_dbh_parents[i];
		for (j in 2:nbYearsPerIndiv[i]) // Loop for all years but the first (which is the parent of indiv i)
		{
			k = k + 1;
			expected_growth = growth(current_dbh, normalised_precip[climate_index[i] + j - 2],
				potentialGrowth, dbh_slope, pr_slope, pr_slope2);
			current_dbh += expected_growth;
			latent_dbh_children[i] += latent_growth[k];
			target += gamma_lpdf(latent_growth[k] | expected_growth^2/processError, expected_growth/processError);
			// j - 2 comes from (j - 1) - 1. The first '-1' is to compensate that j starts at 2, the second '-1' is to match the latent_dbh
		}
	}
	
	// Prior on initial hidden state: This is a diffuse initialisation
	target += uniform_lpdf(latent_dbh_parents | 0.001, 10); // Remember that the dbh is standardised
	
	// Process model
	// for (i in 1:(n_hiddenState - n_indiv))
	// 	target += gamma_lpdf(latent_growth[i] | expected_growth[i]^2/processError, expected_growth[i]/processError);
	
	// --- Observation model
	// Compare true (hidden/latent) parents with observed parents
	target += normal_lpdf(normalised_Yobs[parents_index] | latent_dbh_parents, measureError);

	// Compare true (hidden/latent) children with observed children
	target += normal_lpdf(normalised_Yobs[children_index] | latent_dbh_children, measureError);
}
