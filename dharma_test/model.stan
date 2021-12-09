
data {
	int<lower = 1> n_data;
	int<lower = 1> n_group;
	int<lower = 1> data_per_group[n_group];

	vector[n_data] x;
	vector[n_data] y;
}

parameters {
	// Intercepts
	vector[n_group] beta0_group;
	real beta0;

	// Slope
	real beta1;

	// Error term
	real<lower = 0.00001> sigma_beta;
	real<lower = 0.00001> sigma_res;
}

model {
	// Define variable
	vector[n_data] expected_obs;
	int count = 0;

	// Priors
	target += normal_lpdf(beta0 | 0, 100);
	target += normal_lpdf(beta1 | 0, 100);
	target += normal_lpdf(beta0_group | beta0, sigma_beta);
	target += gamma_lpdf(sigma_beta | 1.0/100, 1.0/100);
	target += gamma_lpdf(sigma_res | 1.0/100, 1.0/100);

	// Model
	for (i in 1:n_group)
	{
		expected_obs[(1 + count):(data_per_group[i] + count)] = beta0_group[i] + beta1*x[(1 + count):(data_per_group[i] + count)];
		count = count + data_per_group[i];
	}

	target += normal_lpdf(y | expected_obs, sigma_res);
}