
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
	array [n_inventories] int<lower = 1, upper = n_children> end_nfi_children; // Ending point of each NFI for children
	
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
	vector [n_plots] plotEffect; // Growth (grouped by plots) when all the explanatory variables are set to 0
	real dbh_slope;

	real pr_slope;
	real pr_slope2;

	real tas_slope;
	real tas_slope2;

	real ph_slope;
	real ph_slope2;

	real competition_slope;

	// Hyper parameters for growth function
	real averageGrowth_mu;
	real<lower = 0> averageGrowth_sd;

	// Errors (observation and process)
	// --- Process error, which is the variance of a gamma distrib /!\
	real<lower = 0.5/sd_dbh^2> sigmaProc;

	// --- Routine observation error, which is constrained by default, see appendix D Eitzel for the calculus.
	array [n_inventories] real<lower = 0.1/sqrt(12)*25.4/sd_dbh> sigmaObs; // Std. Dev. of a normal distrib /!\

	// --- Extreme error, by default at least twice the observation error. Rüger 2011 found it is around 8 times larger
	array [n_inventories] real<lower = 2*0.1/sqrt(12)*25.4/sd_dbh> etaObs; // Std. Dev. of a normal distrib /!\
	
	array [n_inventories] real<lower = 0, upper = 1> proba; // Probabilities of occurrence of extreme errors (etaObs) for each NFI

	// Latent states
	// --- Parent (i.e., primary) dbh
	vector<lower = 0.1/sd_dbh, upper = 3000/sd_dbh>[n_indiv] latent_dbh_parents; // Real (and unobserved) first measurement dbh (parents)

	// --- Growth
	vector<lower = 0>[n_latentGrowth] latent_growth; // Real (and unobserved) yearly growth
}

transformed parameters {
	// Growth (grouped by plots) when all the explanatory variables are set to 0
	vector [n_plots] averageGrowth;
	averageGrowth = averageGrowth_mu + averageGrowth_sd * plotEffect; // <=> averageGrowth ~ normal(averageGrowth_mu, averageGrowth_sd)
}

model {
	// Declare variables
	int growth_counter = 1; // Counter in the latent space
	int children_counter = 1; // Counter in the 'observation space' for chidlren
	int record_children_counter = 1; // Counter to know when we should register the current dbh value into latent_dbh_children
	real expected_growth;
	real temporary;
	vector [n_children] latent_dbh_children;

	// Priors
	// --- Growth parameters
	target += std_normal_lpdf(plotEffect); // <=> normal_lpdf(plotEffect | 0, 1);
	
	target += normal_lpdf(dbh_slope | 0, 5);
	
	target += normal_lpdf(pr_slope | 0, 5);
	target += normal_lpdf(pr_slope2 | 0, 5);

	target += normal_lpdf(tas_slope | 0, 5);
	target += normal_lpdf(tas_slope2 | 0, 5);

	target += normal_lpdf(ph_slope | 0, 5);
	target += normal_lpdf(ph_slope2 | 0, 5);

	target += normal_lpdf(competition_slope | 0, 5);

	// --- Hyper parameters
	target += normal_lpdf(averageGrowth_mu | -4, 10); // sd_dbh * exp(-4) is around 2 to 3 mm, which is a reasonable average growth
	target += gamma_lpdf(averageGrowth_sd | 1.0/100, 1.0/100);

	// --- Errors
	/*
		A note on the two folowwing priors:
		Let m be the mean of the gamma distribution, and v its variance (see comment at the beginning of this file to get shape and rate from
		mean and variance).

		For instance, for the routine measure error, which is the standard deviation of the observation around the latent state, we have:
		m = sqrt(3)
		v = 0.025

		Then, because I work on standardised dbh, I need to divide m by sd(dbh), and v by sd(dbh)^2. Therefore, we get the following:
		shape = (m/sd(dbh))^2 / (v/sd(dbh)^2) = m^2/v
		rate =  (m/sd(dbh))   / (v/sd(dbh)^2) = sd(dbh) * m/v

		For the process error, which is the variance of growth, we have:
		shape = (m/var(dbh))^2 / (v/var(dbh)^2) = m^2/v
		rate = (m/var(dbh))    / (v/var(dbh)^2) = var(dbh) * m/v
	*/
	target += gamma_lpdf(sigmaProc | 5.0^2/1, sd_dbh^2*5.0/1); // Remember that sigmaProc is a variance, not a sd!
	target += gamma_lpdf(sigmaObs | 3.0/0.025, sd_dbh*sqrt(3)/0.025); // <=> routine measurement error (sd) = sqrt(3) mm
	target += gamma_lpdf(etaObs | 25.6^2/6.2, sd_dbh*25.6/6.2); // <=> extreme measurement error (sd) = 25.6 mm

	// Previous prior, from Rüger et al 2011: beta_lpdf(proba | 48.67, 1714.84)
	target += beta_lpdf(proba | 3.95, 391.05); // This corresponds to a 1 % chance extrem error, ± 0.5 %

	// Model
	for (i in 1:n_indiv) // Loop over all the individuals
	{
		temporary = latent_dbh_parents[i];
		for (j in 1:nbYearsGrowth[i]) // Loop for all years but the first (which is the parent of indiv i)
		{
			// Process model
			expected_growth = growth(temporary, normalised_precip[climate_index[i] + j - 1],
				normalised_tas[climate_index[i] + j - 1], normalised_ph[plot_index[i]], normalised_standBasalArea[climate_index[i] + j - 1],
				averageGrowth[plot_index[i]], dbh_slope, pr_slope, pr_slope2, tas_slope, tas_slope2, ph_slope, ph_slope2, competition_slope);
			target += gamma_lpdf(latent_growth[growth_counter] | expected_growth^2/sigmaProc, expected_growth/sigmaProc);

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
		// Compare true (i.e., hidden or latent) parents with observed parents per chunk (i.e. per country/NFI)
		target += (1 - proba[k]) * normal_lpdf(normalised_Yobs[parents_index[start_nfi_parents[k]:end_nfi_parents[k]]] |
				latent_dbh_parents[start_nfi_parents[k]:end_nfi_parents[k]], sigmaObs[k]) +
			proba[k] * normal_lpdf(normalised_Yobs[parents_index[start_nfi_parents[k]:end_nfi_parents[k]]] |
				latent_dbh_parents[start_nfi_parents[k]:end_nfi_parents[k]], etaObs[k]);

		// Compare true (i.e., hidden or latent) children with observed children per chunk (i.e. per country/NFI)
		target += (1 - proba[k]) * normal_lpdf(normalised_Yobs[children_index[start_nfi_children[k]:end_nfi_children[k]]] |
				latent_dbh_children[start_nfi_children[k]:end_nfi_children[k]], sigmaObs[k]) +
			proba[k] * normal_lpdf(normalised_Yobs[children_index[start_nfi_children[k]:end_nfi_children[k]]] |
				latent_dbh_children[start_nfi_children[k]:end_nfi_children[k]], etaObs[k]);
	}
}
