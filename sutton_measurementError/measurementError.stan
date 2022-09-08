
/*
	This model is from Rüger et al. 2011, Growth Strategies of Tropical Tree Species: Disentangling Light and Size Effects
*/

data {
	// Number of data
	int<lower = 1> n_trees; // Number of trees (individuals), also the number of data since each tree is measured exactly twice

	// Data
	array[n_trees] int<lower = 0> dbh1; // dbh measured by person 1 in 2016
	array[n_trees] int<lower = 0> dbh2; // dbh measured by person 2 in 2016
}

parameters {
	// Parameter error
	// --- Routine obs error
	real<lower = 0> sd_a;
	real<lower = 0> sd_b;

	// --- Extreme obs error
	real<lower = 0> etaObs;

	// --- Proba occurrence extreme error
	real<lower = 0, upper = 1> proba;

	// Latent states (dbh)
	array[n_trees] real<lower = 0> latent_dbh;
}

transformed parameters {
	array[n_trees] real sigmaObs;
	for (i in 1:n_trees)
		sigmaObs[i] = sd_a + sd_b*latent_dbh[i];
}

model {
	// Priors
	// --- Routine obs error
	target += gamma_lpdf(sd_a | 0.927^2/0.2, 0.927/0.2); // average from Ruger2011, with 35 times her variance
	target += gamma_lpdf(sd_b | 0.0038^2/4.536e-06, 0.0038/4.536e-06); // average from Ruger2011, with 35 times her variance

	// --- Extreme obs error
	target += gamma_lpdf(etaObs | 30.0^2/225.0, 30.0/225.0); // <=> extreme measurement error (sd) = 30 mm ± 15 mm

	// --- Proba occurrence extreme error
	target += beta_lpdf(proba | 1.453871, 51.22261); // This corresponds to a 2.76 % chance extrem error, ± 2.23 %. Mean from Rüger 2011

	// Diffuse initialisation
	target += uniform_lpdf(latent_dbh | 0, 1000);

	// Likelihood (mixture of distribution)
	for (i in 1:n_trees)
	{
		target += log_mix(proba,
			normal_lpdf(dbh1[i] | latent_dbh[i], etaObs),
			normal_lpdf(dbh1[i] | latent_dbh[i], sigmaObs[i]));

		target += log_mix(proba,
			normal_lpdf(dbh2[i] | latent_dbh[i], etaObs),
			normal_lpdf(dbh2[i] | latent_dbh[i], sigmaObs[i]));
	}
}
