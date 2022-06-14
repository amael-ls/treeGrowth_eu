
/*
	Comments:
		- Yobs is an vector, containing all dbh observations
		- I call "parent" the first measured dbh of each individual. The parents have to be treated separately
		- The observation error for the French data is different from the other countries because there have been corrections in the data.
			Therefore, there is no extrem error due to, for instance, typos.
		- Note that the dbh is standardised but not centred. The reason was that in a previous version of the model, it would have
			implied a division by 0. Since then, I do not have this problem anymore, but I kept the dbh non-centred.
	
	Reminder (method of first and second moments, i.e., mean and variance or standard deviation):
		- The gamma distribution of stan uses, in this order, shape (alpha) and rate (beta). It can however be reparametrised
			using mean and var: alpha = mean^2/var, beta = mean/var
		
		- The lognormal distribution of stan uses, in this order, meanlog (mu) and sdlog (sigma). It can however be reparametrised
			using mean and sd: mu = log(mean^2/sqrt(sd^2 + mean^2)), sigma = sqrt(log(sd^2/mean^2 + 1))
*/

functions {
	// Function for growth. This returns the expected growth in mm, for 1 year.
	real growth(real dbh0, real precip, real temperature, real ph, real standBasalArea, real averageGrowth, real dbh_slope,
		real pr_slope, real pr_slope2 , real tas_slope, real tas_slope2, real ph_slope, real ph_slope2, real competition_slope)
	{
		/*
			It takes the following variables:
				- dbh0: The diameter
				- precip: The precipitations
				- temperature: The temperature
				- ph: The pH (CaCl2)
				- standBasalArea: Basal area of the stand (i.e., includes all the trees, regardless of the species)
			
			It takes the following parameters:
				- averageGrowth: The basic growth, when all the explanatory variables are set to zero
				- dbh_slope: The slope for dbh
				- pr_slope: The slope for precipitation
				- pr_slope2: The slope for precipitation (squared term)
				- tas_slope: The slope for temperature (tas)
				- tas_slope2: The slope for temperature (tas, squared term)
				- ph_slope: The slope for pH
				- ph_slope2: The slope for pH (squared term)
				- competition_slope: The slope for competition
		*/

		return (exp(averageGrowth + dbh_slope*dbh0 + pr_slope*precip + pr_slope2*precip^2 +
			tas_slope*temperature + tas_slope2*temperature^2 + ph_slope*ph + ph_slope2*ph^2 + competition_slope*standBasalArea));
	}
}

data {
	// Number of data
	int<lower = 1> n_indiv; // Total number of individuals (all NFIs together)
	int<lower = 1> n_climate; // Dimension of the climate vector (all NFIs together)
	int<lower = 1, upper = n_indiv> n_plots; // Number of plots (all NFIs together)
	int<lower = 1> n_obs; // Total number of tree observations (all NFIs together)
	int<lower = 1> n_latentGrowth; // Dimension of the state space vector for latent growth (all NFIs together)
	int<lower = n_obs - n_indiv, upper = n_obs - n_indiv> n_children; // Number of children tree observations = n_obs - n_indiv
	array [n_indiv] int<lower = 2, upper = n_obs> nbYearsGrowth; // Number of years of growth for each individual
	int<lower = 1> n_inventories; // Number of forest inventories involving different measurement errors in the data

	// Indices
	array [n_indiv] int<lower = 1, upper = n_obs - 1> parents_index; // Index of each parent in the 'observation space'
	array [n_children] int<lower = 2, upper = n_obs> children_index; // Index of children in the 'observation space'
	array [n_children] int<lower = 2> latent_children_index; // Index of children in the 'latent space'
	array [n_indiv] int<lower = 1, upper = n_climate - 1> climate_index; // Index of the climate associated to each parent
	
	array [n_inventories] int<lower = 1, upper = n_indiv - 1> start_nfi_parents; // Starting point of each NFI for parents
	array [n_inventories] int<lower = 2, upper = n_indiv> end_nfi_parents; // Ending point of each NFI for parents
	array [n_inventories] int<lower = 1, upper = n_children - 1> start_nfi_children; // Starting point of each NFI for children
	array [n_inventories] int<lower = 2, upper = n_children> end_nfi_children; // Ending point of each NFI for children
	
	array [n_indiv] int<lower = 1, upper = n_plots> plot_index; // Indicates to which plot individuals belong to

	// Observations
	vector<lower = 0>[n_obs] Yobs;

	// sd_dbh is for the whole species and should not be more than 5% different from the sd of the subsample, namely sd(Yobs)
	real<lower = 0.95*sd(Yobs), upper = 1.05*sd(Yobs)> sd_dbh;

	// Explanatory variables
	vector<lower = 0>[n_climate] precip; // Precipitations
	real<lower = 0> pr_mu; // To standardise the precipitations
	real<lower = 0> pr_sd; // To standardise the precipitations

	vector[n_climate] tas; // Temperature
	real tas_mu; // To standardise the temperature
	real<lower = 0> tas_sd; // To standardise the temperature

	vector<lower = 0, upper = 14>[n_plots] ph; // pH of the soil measured with CaCl2
	real<lower = 0, upper = 14> ph_mu; // To standardise the pH
	real<lower = 0> ph_sd; // To standardise the pH

	vector<lower = 0>[n_climate] standBasalArea; // Sum of the tree basal area for a given plot at a given time (interpolation for NA data)
	real ba_mu; // To standardise the basal area
	real<lower = 0> ba_sd; // To standardise the basal area
}

transformed data {
	vector[n_obs] normalised_Yobs = Yobs/sd_dbh; // Normalised but NOT centred dbh
	vector[n_climate] normalised_precip = (precip - pr_mu)/pr_sd; // Normalised and centred precipitations
	vector[n_climate] normalised_tas = (tas - tas_mu)/tas_sd; // Normalised and centred temperatures
	vector[n_plots] normalised_ph = (ph - ph_mu)/ph_sd; // Normalised and centred pH
	vector[n_climate] normalised_standBasalArea = (standBasalArea - ba_mu)/ba_sd; // Normalised and centred BA
}

parameters {
	// Parameters for growth function
	real averageGrowth;
	real dbh_slope;

	real pr_slope;
	real pr_slope2;

	real tas_slope;
	real tas_slope2;

	real ph_slope;
	real ph_slope2;

	real competition_slope;
	
	// Errors (observation and process)
	// --- Process error, which is the sd of a lognormal distrib /!\
	// real<lower = 0.5/sd_dbh^2> sigmaProc;

	// --- Routine observation error, which is constrained by default, see appendix D Eitzel for the calculus.
	// array [n_inventories] real<lower = 0.1/sqrt(12)*25.4/sd_dbh> sigmaObs; // Std. Dev. of a normal distrib /!\

	// --- Extreme error, by default at least twice the observation error. Rüger 2011 found it is around 8 times larger
	// array [n_inventories] real<lower = 2*0.1/sqrt(12)*25.4/sd_dbh> etaObs; // Std. Dev. of a normal distrib /!\
	
	array [n_inventories] real<lower = 0, upper = 1> proba; // Probabilities of occurrence of extreme errors (etaObs) for each NFI

	// Latent states
	// --- Parent (i.e., primary) dbh
	vector<lower = 0.1/sd_dbh, upper = 3000/sd_dbh>[n_indiv] latent_dbh_parents; // Real (and unobserved) first measurement dbh (parents)

	// --- Growth
	vector<lower = 0>[n_latentGrowth] latent_growth; // Real (and unobserved) yearly growth
}

model {
	// Declare variables
	int growth_counter = 1; // Counter in the latent space
	int children_counter = 1; // Counter in the 'observation space' for chidlren
	int record_children_counter = 1; // Counter to know when we should register the current dbh value into latent_dbh_children
	real expected_growth;
	real temporary;
	vector [n_children] latent_dbh_children;
	real growth_mean_logNormal;
	real growth_sd_logNormal;

	real sigmaProc = 0.008977944;
	real sigmaObs = 0.01209929;
	real etaObs = 0.1788352;

	// Priors
	// --- Growth parameters
	target += normal_lpdf(averageGrowth | -4, 10); // sd_dbh * exp(-4) is around 2 to 3 mm (reasonable average growth) for sd_dbh = 130
	target += normal_lpdf(dbh_slope | 0, 5);
	
	target += normal_lpdf(pr_slope | 0, 5);
	target += normal_lpdf(pr_slope2 | 0, 5);

	target += normal_lpdf(tas_slope | 0, 5);
	target += normal_lpdf(tas_slope2 | 0, 5);

	target += normal_lpdf(ph_slope | 0, 5);
	target += normal_lpdf(ph_slope2 | 0, 5);

	target += normal_lpdf(competition_slope | 0, 5);

	// --- Errors
	/*
		A note on the folowing priors:
		--> Let m be the mean of the gamma distribution, and v its variance (see comment at the beginning of this file to get shape and
			rate from mean and variance).

		For instance, for the routine measure error, which is the standard deviation of the observation around the latent state, we have:
		m = sqrt(3)
		v = 0.025

		Then, because I work on standardised dbh, I need to divide m by sd(dbh), and v by sd(dbh)^2. Therefore, we get the following:
		shape = (m/sd(dbh))^2 / (v/sd(dbh)^2) = m^2/v
		rate =  (m/sd(dbh))   / (v/sd(dbh)^2) = sd(dbh) * m/v

		--> Let m be the mean of the lognormal distribution, and s its standard deviation (see comment at the beginning of this file to get
		meanlog and sdlog from mean and sd).
		For the process error, which is the sd of growth, we have:
		meanlog = log( (m/sd(dbh))^2/sqrt((s/sd(dbh))^2 + (m/sd(dbh))^2) ) = log(m^2/sqrt(s^2 + m^2)) - log(sd(dbh)) = meanlog - log(sd(dbh))
		sdlog = sqrt(log( (s/sd(dbh))^2 / (m/sd(dbh))^2 + 1)) = sqrt(log(s^2/m^2 + 1))

		However, for the lognormal, I used meanlog = log(mu_h) = log(1.28), and sd_log = sigma_h = 0.16 (Nadja Rüger, personnal
		communication)

		The values are taken from Rüger et al (2011), Growth Strategies of Tropical Tree Species: Disentangling Light and Size Effects
		for sigmaProc (personnal communication), etaObs and proba priors.

		The values are taken from Luoma et al (2017), Assessing Precision in Conventional Field Measurements of Individual Tree Attributes
		for sigmaObs prior.

		The average values of the priors on the standardised scale are approximately (obtained with one sample of 1e8):
		sigmaProc: 0.008977944
		sigmaObs: 0.01209929
		etaObs: 0.1788352
	*/
	// target += lognormal_lpdf(sigmaProc | 0.2468601 - log(sd_dbh), 0.09); // <=> procError = 1.30 mm/yr ± 0.21 mm/yr
	// target += gamma_lpdf(sigmaObs | 3.0/0.025, sd_dbh*sqrt(3)/0.025); // <=> routine measurement error (sd) = sqrt(3) mm
	// target += gamma_lpdf(etaObs | 25.6^2/6.2, sd_dbh*25.6/6.2); // <=> extreme measurement error (sd) = 25.6 mm

	// target += beta_lpdf(proba | 3.95, 391.05); // This corresponds to a 1 % chance extrem error, ± 0.5 %
	target += beta_lpdf(proba | 48.67, 1714.84); // This corresponds to a 2.76 % chance extrem error, ± 0.39 % (Rüger et. al 2011)

	// Model
	for (i in 1:n_indiv) // Loop over all the individuals
	{
		temporary = latent_dbh_parents[i];
		for (j in 1:nbYearsGrowth[i]) // Loop for all years but the first (which is the parent of indiv i)
		{
			// Process model
			expected_growth = growth(temporary, normalised_precip[climate_index[i] + j - 1],
				normalised_tas[climate_index[i] + j - 1], normalised_ph[plot_index[i]], normalised_standBasalArea[climate_index[i] + j - 1],
				averageGrowth, dbh_slope, pr_slope, pr_slope2, tas_slope, tas_slope2, ph_slope, ph_slope2, competition_slope);

			growth_mean_logNormal = log(expected_growth^2/sqrt(sigmaProc^2 + expected_growth^2));
			growth_sd_logNormal = sqrt(log(sigmaProc^2/expected_growth^2 + 1));

			target += lognormal_lpdf(latent_growth[growth_counter] | growth_mean_logNormal, growth_sd_logNormal);

			// Dbh at time t + 1
			temporary += latent_growth[growth_counter];

			growth_counter += 1;
			record_children_counter += 1;

			// Only the relevant (i.e., children) dbh are recorded
			if (record_children_counter == latent_children_index[children_counter])
			{
				latent_dbh_children[children_counter] = temporary;
				children_counter += 1;
			}
		}
		record_children_counter += 1; // The growth counter stops one year earlier! Indeed 3 measurements implies only 2 growing years!
	}
	
	// Prior on initial hidden state: This is a diffuse initialisation
	target += uniform_lpdf(latent_dbh_parents | 0.1/sd_dbh, 3000/sd_dbh); // Remember that the dbh is in mm and standardised
	
	// --- Observation model
	for (k in 1:n_inventories)
	{
		// Compare true (i.e., hidden or latent) parents with observed parents
		// Do not try to vectorise here! https://mc-stan.org/docs/2_29/stan-users-guide/vectorizing-mixtures.html
		for (i in start_nfi_parents[k]:end_nfi_parents[k])
			target += log_mix(proba[k], normal_lpdf(normalised_Yobs[parents_index[i]] | latent_dbh_parents[i], etaObs),
				normal_lpdf(normalised_Yobs[parents_index[i]] | latent_dbh_parents[i], sigmaObs));

		// Compare true (i.e., hidden or latent) children with observed children
		for (i in start_nfi_children[k]:end_nfi_children[k])
			target += log_mix(proba[k], normal_lpdf(normalised_Yobs[children_index[i]] | latent_dbh_children[i], etaObs),
				normal_lpdf(normalised_Yobs[children_index[i]] | latent_dbh_children[i], sigmaObs));
	}
}
