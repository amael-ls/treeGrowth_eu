
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
	real growth(real dbh0, real precip, real temperature, real standBasalArea,
		real potentialGrowth, real dbh_slope, real pr_slope, real pr_slope2 , real tas_slope, real tas_slope2, real competition_slope)
	{
		/*
			It takes the following variables:
				- dbh0: The diameter
				- precip: The precipitations
				- temperature: The temperature
				- standBasalArea: Basal area of the stand (i.e., includes all the trees, regardless of the species)
			
			It takes the following parameters:
				- potentialGrowth: The basic growth, when all the explanatory variables are set to zero
				- dbh_slope: The slope for dbh
				- pr_slope: The slope for precipitation
				- pr_slope2: The slope for precipitation (squared term)
				- tas_slope: The slope for temperature (tas)
				- tas_slope2: The slope for temperature (tas, squared term)
				- competition_slope: The slope for competition
		*/

		return (exp(potentialGrowth + dbh_slope*dbh0 + pr_slope*precip + pr_slope2*precip^2 +
			tas_slope*temperature + tas_slope2*temperature^2 + competition_slope*standBasalArea));
	}
}

data {
	// Number of data
	int<lower = 1> n_indiv; // Total number of individuals (all NFIs together)
	int<lower = 1> n_climate; // Dimension of the climate vector (all NFIs together)
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

	// Observations
	vector<lower = 0>[n_obs] Yobs;

	// sd_dbh is for the whole species and should not be more than 5% different from the sd of the subsample, namely sd(Yobs)
	real<lower = 0.95*sd(Yobs), upper = 1.05*sd(Yobs)> sd_dbh;

	// Explanatory variables
	vector<lower = 0>[n_climate] precip; // Precipitations
	real pr_mu; // To standardise the precipitations
	real<lower = 0> pr_sd; // To standardise the precipitations

	vector[n_climate] tas; // Temperature
	real tas_mu; // To standardise the temperature
	real<lower = 0> tas_sd; // To standardise the temperature

	/*
		The tree weight correspond to the weight of the tree per hectar. It accounts for:
			1. the area of the plot (disc with a radius of 6, 9 or 15 metres depending on the size class)
			2. the proximity to the boundary (rectification of the probability of drawing trees at the edge)
	*/
	vector<lower = 0>[n_indiv] standBasalArea; // Sum of the tree basal area for a given plot at a given time
}

transformed data {
	vector[n_obs] normalised_Yobs = Yobs/sd_dbh; // Normalised but NOT centred dbh
	vector[n_climate] normalised_precip = (precip - pr_mu)/pr_sd; // Normalised and centred precipitations
	vector[n_climate] normalised_tas = (tas - tas_mu)/tas_sd; // Normalised and centred temperatures
	vector[n_indiv] normalised_standBasalArea = (standBasalArea - mean(standBasalArea))/sd(standBasalArea); // Normalised and centred BA
}

parameters {
	// Parameters
	real potentialGrowth; // Growth when all the explanatory variables are set to 0
	real dbh_slope;

	real pr_slope;
	real pr_slope2;

	real tas_slope;
	real tas_slope2;

	real competition_slope;

	real<lower = 0.5/sd_dbh^2> sigmaProc; // /!\ Variance of a gamma distrib /!\
	// The observation error is constrained by default, see appendix D Eitzel for the calculus.
	array [n_inventories] real<lower = 0.1/sqrt(12)*25.4/sd_dbh> sigmaObs; // /!\ Std. Dev. of a normal distrib /!\
	// I decided that the extreme error would be at least twice the observation error. Rüger 2011 found it is around 8 times larger
	array [n_inventories] real<lower = 2*0.1/sqrt(12)*25.4/sd_dbh> etaObs; // /!\ Std. Dev. of a normal distrib /!\
	
	array [n_inventories] real<lower = 0, upper = 1> p; // Probabilities of occurrence of extreme errors (etaObs) for each NFI

	vector<lower = 0.1/sd_dbh, upper = 3000/sd_dbh>[n_indiv] latent_dbh_parents; // Real (and unobserved) first measurement dbh (parents)
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

	// Priors
	target += normal_lpdf(potentialGrowth | 0, 100);
	target += normal_lpdf(dbh_slope | 0, 5);
	
	target += normal_lpdf(pr_slope | 0, 5);
	target += normal_lpdf(pr_slope2 | 0, 5);

	target += normal_lpdf(tas_slope | 0, 5);
	target += normal_lpdf(tas_slope2 | 0, 5);

	target += normal_lpdf(competition_slope | 0, 5);

	/*
		A note on the two folowwing priors:
		Let m be the mean of the gamma distribution, and v its variance (see comment at the beginning of this file to get shape and rate from
		mean and variance).

		For instance, for the routine measure error, which is the standard deviation of the observation around the latent state, we have:
		m = sqrt(3)
		v = 0.15

		Then, because I work on standardised dbh, I need to divide m by sd(dbh), and v by sd(dbh)^2. Therefore, we get the following:
		shape = (m/sd(dbh))^2 / (v/sd(dbh)^2) = m^2/v
		rate =  (m/sd(dbh))   / (v/sd(dbh)^2) = sd(dbh) * m/v

		For the process error, which is the variance of growth, we have:
		shape = (m/var(dbh))^2 / (v/var(dbh)^2) = m^2/v
		rate = (m/var(dbh))    / (v/var(dbh)^2) = var(dbh) * m/v
	*/
	target += gamma_lpdf(sigmaProc | 5.0^2/1, variance(Yobs)*5.0/1); // Remember that sigmaProc is a variance, not a sd!
	target += gamma_lpdf(sigmaObs | 3.0/0.15, sd_dbh*sqrt(3)/0.15); // <=> routine measurement error (sd) = sqrt(3) mm
	target += gamma_lpdf(etaObs | 25.6^2/6.2, sd_dbh*25.6/6.2); // <=> extreme measurement error (sd) = 25.6 mm

	target += beta_lpdf(p | 48.67, 1714.84);

	// Model
	for (i in 1:n_indiv) // Loop over all the individuals
	{
		temporary = latent_dbh_parents[i];
		for (j in 1:nbYearsGrowth[i]) // Loop for all years but the first (which is the parent of indiv i)
		{
			// Process model
			expected_growth = growth(temporary, normalised_precip[climate_index[i] + j - 1],
				normalised_tas[climate_index[i] + j - 1], normalised_standBasalArea[i],
				potentialGrowth, dbh_slope, pr_slope, pr_slope2, tas_slope, tas_slope2, competition_slope);
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
		target += (1 - p[k]) * normal_lpdf(normalised_Yobs[parents_index[start_nfi_parents[k]:end_nfi_parents[k]]] |
				latent_dbh_parents[start_nfi_parents[k]:end_nfi_parents[k]], sigmaObs[k]) +
			p[k] * normal_lpdf(normalised_Yobs[parents_index[start_nfi_parents[k]:end_nfi_parents[k]]] |
				latent_dbh_parents[start_nfi_parents[k]:end_nfi_parents[k]], etaObs[k]);

		// Compare true (i.e., hidden or latent) children with observed children per chunk (i.e. per country/NFI)
		target += (1 - p[k]) * normal_lpdf(normalised_Yobs[children_index[start_nfi_children[k]:end_nfi_children[k]]] |
				latent_dbh_children[start_nfi_children[k]:end_nfi_children[k]], sigmaObs[k]) +
			p[k] * normal_lpdf(normalised_Yobs[children_index[start_nfi_children[k]:end_nfi_children[k]]] |
				latent_dbh_children[start_nfi_children[k]:end_nfi_children[k]], etaObs[k]);
	}
}
