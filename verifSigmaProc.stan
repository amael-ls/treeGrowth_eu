
functions {
	real growth_fct(real dbh0, real averageGrowth, real dbh_slope)
	{
		return (exp(averageGrowth + dbh_slope*dbh0));
	}
}

data {
	int<lower = 0> n;
	int<lower = 0> n_indiv;
	real<lower = 0> sd_dbh;

	array [n_indiv] int parents_index;
	array [n_indiv] int children_index;

	vector[n] dbh;
}

transformed data {
	vector[n] normalised_dbh = dbh/sd_dbh; // Normalised but NOT centred dbh
	vector [n_indiv] growth = normalised_dbh[children_index] - normalised_dbh[parents_index];
}

parameters {
	real averageGrowth;
	real dbh_slope;

	real<lower = 0> sigmaProc;
}

model {
	real expected_growth;
	real growth_mean_logNormal;
	real growth_sd_logNormal;
	
	target += normal_lpdf(averageGrowth | 0, 5);
	target += normal_lpdf(dbh_slope | 0, 5);

	target += lognormal_lpdf(sigmaProc | 1.256145 - log(sd_dbh), 0.2696117); // <=> procError = 3.64 mm/yr ± 1 mm/yr

	for (i in 1:n_indiv)
	{
		expected_growth = growth_fct(normalised_dbh[parents_index[i]], averageGrowth, dbh_slope);
	
		growth_mean_logNormal = log(expected_growth^2/sqrt(sigmaProc^2 + expected_growth^2));
		growth_sd_logNormal = sqrt(log(sigmaProc^2/expected_growth^2 + 1));

		target += lognormal_lpdf(growth[i] | growth_mean_logNormal, growth_sd_logNormal);
	}
}
