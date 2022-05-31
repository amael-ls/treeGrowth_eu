
/*
	Comments:
	Script to generate new data. Note that 'model' is an unnecessary block here as restart from results. gq stands for generated quantities
	More information can be found at https://mc-stan.org/cmdstanr/reference/model-method-generate-quantities.html
	
	The 'new observations' (poserior predictions) are simulated as follow:
		1. Draw the vector of parameters theta (which includes the latent states!)
		2. Generate the parent observation from the corresponding latent dbh according to the model
		3. Generate the children observation according to the model (i.e., after n years of yearly latent growth)
*/

functions {
	// Function for growth. This returns the expected growth in mm, for 1 year.
	real growth(real dbh0, real precip, real temperature, real ph, real standBasalArea, real averageGrowth, real dbh_slope,
		real pr_slope, real pr_slope2 , real tas_slope, real tas_slope2, real ph_slope, real ph_slope2, real competition_slope)
	{
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
	
	array [n_indiv] int<lower = 1, upper = n_inventories> nfi_id; // Indicates to which NFI individuals belong to
	
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

generated quantities {
	// Declaration output variables
	array [n_obs] real newObservations;
	array [n_obs] real latent_dbh_parentsChildren;
	array [n_obs] real proba_extreme_obs; // Probability that obsertation i is from the extreme part of the mixture
	array [n_latentGrowth] real latentG_residuals;
	array [n_latentGrowth + n_indiv] real yearly_latent_dbh;
	simplex [2] proba_small_large;
	
	{
		// Variables declared in nested blocks are local variables, not generated quantities, and thus won't be printed.
		real current_latent_dbh;
		int growth_counter = 1;
		int dbh_counter = 1;
		int children_counter = 1;
		int groupError; // Random variable to select from which distribution in the mixture the data is from
		array [2] real obsError_small_large;
		real expected_growth;
		array [4]  real shitty_array;

		array [n_children] real temporary_children; // To store the children, will be added to latent_dbh_parentsChildren
		array [n_children] real temporary_observation; // To store the children, will be added to newObservations

		for (i in 1:n_indiv) // Loop over all the individuals
		{
			// Fill the parents
			latent_dbh_parentsChildren[parents_index[i]] = latent_dbh_parents[i];
			yearly_latent_dbh[dbh_counter] = latent_dbh_parents[i];

			// Generate the parent observation conditional on the parent state
			proba_small_large[1] = 1 - proba[nfi_id[i]]; // Remember that proba is the probability of committing a LARGE error
			proba_small_large[2] = proba[nfi_id[i]];
			obsError_small_large[1] = sigmaObs[nfi_id[i]];
			obsError_small_large[2] = etaObs[nfi_id[i]];
			groupError = categorical_rng(proba_small_large);
			newObservations[parents_index[i]] = normal_rng(latent_dbh_parents[i], obsError_small_large[groupError]);

			// Compute the probability that the parent observation is from extreme part of the mixture
			shitty_array[1] = normal_lpdf(normalised_Yobs[parents_index[i]] | latent_dbh_parents[i], obsError_small_large[1]);
			shitty_array[2] = normal_lpdf(normalised_Yobs[parents_index[i]] | latent_dbh_parents[i], obsError_small_large[2]);
			shitty_array[3] = gamma_lpdf(obsError_small_large[1] | 3.0/0.025, sd_dbh*sqrt(3)/0.025); //! WRONG, ...
			//! ... should be [s_n = 1 | obsError_small_large[1]], check
			//! https://mc-stan.org/docs/2_29/functions-reference/categorical-distribution.html
			//! There is also softmax, check what it is
			shitty_array[4] = gamma_lpdf(obsError_small_large[2] | 25.6^2/6.2, sd_dbh*25.6/6.2); //! WRONG
			proba_extreme_obs[parents_index[i]] =
				normal_lpdf(normalised_Yobs[parents_index[i]] | latent_dbh_parents[i], obsError_small_large[2]) +
				gamma_lpdf(obsError_small_large[2] | 25.6^2/6.2, sd_dbh*25.6/6.2) -  //! WRONG
				log_sum_exp(shitty_array);

			// Starting point to compute latent dbh child
			current_latent_dbh = latent_dbh_parents[i];

			for (j in 1:nbYearsGrowth[i]) // Loop for all years but the first (which is the parent of indiv i)
			{
				// Process model
				expected_growth = growth(current_latent_dbh, normalised_precip[climate_index[i] + j - 1],
					normalised_tas[climate_index[i] + j - 1], normalised_ph[plot_index[i]],
					normalised_standBasalArea[climate_index[i] + j - 1], averageGrowth[plot_index[i]],
					dbh_slope, pr_slope, pr_slope2, tas_slope, tas_slope2, ph_slope, ph_slope2, competition_slope);

				// Record difference expected growth minus latent growth
				latentG_residuals[growth_counter] = expected_growth - latent_growth[growth_counter];

				// Dbh at time t + 1
				current_latent_dbh += latent_growth[growth_counter]; // Or should it be += gamma(mean = latent_growth, var = sigmaProc)?
				growth_counter += 1;
				dbh_counter += 1;
				yearly_latent_dbh[dbh_counter] = current_latent_dbh;

				// Recording the children dbh
				if (dbh_counter == latent_children_index[children_counter])
				{
					temporary_children[children_counter] = current_latent_dbh;
					groupError = categorical_rng(proba_small_large);
					temporary_observation[children_counter] = normal_rng(current_latent_dbh, obsError_small_large[groupError]);

					shitty_array[1] =
						normal_lpdf(normalised_Yobs[children_index[children_counter]] | current_latent_dbh, obsError_small_large[1]);
					shitty_array[2] =
						normal_lpdf(normalised_Yobs[children_index[children_counter]] | current_latent_dbh, obsError_small_large[2]);
					shitty_array[3] = gamma_lpdf(obsError_small_large[1] | 3.0/0.025, sd_dbh*sqrt(3)/0.025); //! WRONG
					shitty_array[4] = gamma_lpdf(obsError_small_large[2] | 25.6^2/6.2, sd_dbh*25.6/6.2); //! WRONG

					proba_extreme_obs[children_index[children_counter]] =
					normal_lpdf(normalised_Yobs[children_index[children_counter]] | current_latent_dbh, obsError_small_large[2]) +
					gamma_lpdf(obsError_small_large[2] | 25.6^2/6.2, sd_dbh*25.6/6.2) -  //! WRONG
					log_sum_exp(shitty_array);
					children_counter += 1;
				}
			}
			dbh_counter += 1;
		}
		// Merging
		latent_dbh_parentsChildren[children_index] = temporary_children;
		newObservations[children_index] = temporary_observation;
	}
}
