
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
	int<lower = 1> n_latentGrowth; // Dimension of the state space vector for latent growth
	int<lower = n_obs - n_indiv, upper = n_obs - n_indiv> n_children; // Number of children trees observations = n_obs - n_indiv
	array [n_indiv] int<lower = 2, upper = n_obs> nbYearsGrowth; // Number of years of growth for each individual

	// Indices
	array [n_indiv] int<lower = 1, upper = n_obs - 1> parents_index; // Index of each parent in the observed data
	array [n_children] int<lower = 2, upper = n_obs> children_index; // Index of children in the observed data
	array [n_indiv] int<lower = 1, upper = n_precip - 1> climate_index; // Index of the climate associated to each parent

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

	real<lower = 0.5/sd(Yobs)^2> processError;
	real<lower = 0.1/sqrt(12)*25.4/sd(Yobs)> measureError; // Constrained by default, see appendix D Eitzel for the calculus

	vector<lower = 0, upper = 10>[n_indiv] latent_dbh_parents; // Real (and unobserved) first measurement dbh (parents)
	vector<lower = 0>[n_latentGrowth] latent_growth; // Real (and unobserved) yearly growth
}

model {
	// Declare variables
	int growth_counter = 1;
	real expected_growth;
	// real processError = 5.68/sd(Yobs)^2; // /!\ This is a variance, contrary to measureError which is a std. dev.
	vector [n_obs - n_indiv] latent_dbh_children;

	// Priors
	target += normal_lpdf(potentialGrowth | 0, 100);
	target += normal_lpdf(dbh_slope | 0, 5);
	
	target += normal_lpdf(pr_slope | 0, 5);
	target += normal_lpdf(pr_slope2 | 0, 5);

	/*
		A note on the two folowwing priors:
		Let m be the mean of the gamma distribution, and v its variance (see comment at the beginning of this file to get shape and rate from
		mean and variance).

		For instance, for the measurement error, which is the standard deviation of the observation around the latent state, we have:
		m = sqrt(3)
		v = 0.15

		Then, because I work on standardised dbh, I need to divide m by sd(dbh), and v by sd(dbh)^2. Therefore, we get the following:
		shape = (m/sd(dbh))^2 / (v/sd(dbh)^2) = m^2/v
		rate =  (m/sd(dbh))   / (v/sd(dbh)^2) = sd(dbh) * m/v

		For the process error, which is the variance of growth, we have:
		shape = (m/var(dbh))^2 / (v/var(dbh)^2) = m^2/v
		rate = (m/var(dbh))    / (v/var(dbh)^2) = var(dbh) * m/v
	*/
	target += gamma_lpdf(processError | 5.0^2/1, variance(Yobs)*5.0/1); // Remember that processError is a variance, not a sd!
	target += gamma_lpdf(measureError | 3.0/0.001, sd(Yobs)*sqrt(3)/0.001); // <=> measurement error (sd) = sqrt(3) mm

	// Model
	for (i in 1:n_indiv) // Loop over all the individuals
	{
		latent_dbh_children[i] = latent_dbh_parents[i];
		for (j in 1:nbYearsGrowth[i]) // Loop for all years but the first (which is the parent of indiv i)
		{
			// Process model
			expected_growth = growth(latent_dbh_children[i], normalised_precip[climate_index[i] + j - 1],
				potentialGrowth, dbh_slope, pr_slope, pr_slope2);
			target += gamma_lpdf(latent_growth[growth_counter] | expected_growth^2/processError, expected_growth/processError);

			// Dbh at time t + 1, only the last (i.e., child) dbh is recorded
			latent_dbh_children[i] += latent_growth[growth_counter];
			growth_counter += 1;
		}
	}
	
	// Prior on initial hidden state: This is a diffuse initialisation
	target += uniform_lpdf(latent_dbh_parents | 0.001, 10); // Remember that the dbh is standardised
	
	// --- Observation model
	// Compare true (i.e., hidden or latent) parents with observed parents
	target += normal_lpdf(normalised_Yobs[parents_index] | latent_dbh_parents, measureError);

	// Compare true (i.e., hidden or latent) children with observed children
	target += normal_lpdf(normalised_Yobs[children_index] | latent_dbh_children, measureError);
}
