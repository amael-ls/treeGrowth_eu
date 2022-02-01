
/*
	I am using the beta binomial distribution. If I have not done any mistake, for n > 1 fixed, the function relating mean and variance to
		the parameters alpha and beta is bijective. That is to say, I can fix arbitrarily the integer n to a value at least 2 and I am still
		sure that I cover the whole space which is 'the upper right square of the plan' (i.e., the product of positive reals).
*/

data {
	// Number of data
	int<lower = 1> n_trees; // Number of trees

	// Data
	int<lower = 0> dbh1[n_trees]; // dbh measured by person 1
	int<lower = 0> dbh2[n_trees]; // dbh measured by person 2
}

parameters {
	// Parameter error
	real<lower = 0> error;

	// Latent state
	vector<lower = 0>[n_trees] latent_dbh;
}

model {
	// Defines variables
	target += gamma_lpdf(error | 3.0^2/1, 3.0/1); // Gives an average of 3 and a variance of 1

	for (i in 1:n_trees)
	{
		target += normal_lpdf(latent_dbh[i] | dbh1[i], 5);
		target += normal_lpdf(dbh1[i] | latent_dbh[i], error);
		target += normal_lpdf(dbh2[i] | latent_dbh[i], error);
	}
}

// generated quantities {
// 	vector[n_trees] latent_dbh;
// 	vector[n_trees] sigma;

// 	for (i in 1:n_trees)
// 	{
// 		latent_dbh[i] = N*x/(alpha[i] + beta[i])
// 		sigma[i] = N*alpha[i]*beta[i]*(alpha[i] + beta[i] + N)/((alpha[i] + beta[i])^2*(alpha[i] + beta[i] + 1))
// 	}
// }


// parameters {
// 	// Parameter error
// 	real<lower = 0.5/sqrt(12)> error; // The error is at least the error of the tool. See appendix D Eitzel for the calculus
	
// 	// Latent states
// 	vector<lower = 0, upper = N>[n_trees] latent_dbh; // The average dbh
// }

// model {
// 	// Defines variables
// 	real alpha;
// 	real beta;

// 	// Diffuse initialisation of the states
// 	target += uniform_lpdf(latent_dbh | 0, N);

// 	// Error prior
// 	target += gamma_lpdf(error | 3.0^2/100, 3.0/100);

// 	for (i in 1:n_trees)
// 	{
// 		alpha = (-latent_dbh[i]^3 + latent_dbh[i]^2*N - latent_dbh[i]*error)/(latent_dbh[i]^2 - latent_dbh[i]*N + error*N);
// 		if (alpha < 0)
// 		{
// 			print("latent = ", latent_dbh[i]);
// 			print("num = ", -latent_dbh[i]^3 + latent_dbh[i]^2*N - latent_dbh[i]*error);
// 			print("denom = ", latent_dbh[i]^2 - latent_dbh[i]*N + error*N);
// 		}
// 		beta = (latent_dbh[i] - N)*(latent_dbh[i]^2 - latent_dbh[i]*N + error)/(latent_dbh[i]^2 - latent_dbh[i]*N + error*N);

// 		target += beta_binomial_lpmf(dbh1[i] | N, alpha, beta);
// 		target += beta_binomial_lpmf(dbh2[i] | N, alpha, beta);
// 	}
// }
