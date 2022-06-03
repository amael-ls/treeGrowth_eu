
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
	vector [n_indiv] expected_growth;
	target += normal_lupdf(averageGrowth | 0, 5);
	target += normal_lupdf(dbh_slope | 0, 5);

	target += gamma_lupdf(sigmaProc | 3.64^2/1.5, sd_dbh^2*3.64/1.5);

	for (i in 1:n_indiv)
		expected_growth[i] = growth_fct(normalised_dbh[parents_index[i]], averageGrowth, dbh_slope);
	
	target += gamma_lupdf(growth | expected_growth^2/sigmaProc, expected_growth/sigmaProc);
}
