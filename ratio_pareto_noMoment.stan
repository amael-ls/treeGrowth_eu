

//// Comments
// Ref: http://www.naturalspublishing.com/files/published/73t3w33ckp06vv.pdf
//	On Three Parameter Weighted Pareto Type II Distribution: Properties and Applications in Medical Sciences
//
// Based on the aforementioned ref, I think it would makes sense that alpha is an increasing function of dbh. This gives:
// 	- A decreasing average in function of dbh (towards 1 if distribution shifted, as it is the case in Pareto II)
//	- A decreasing variance in function of dbh
// The other parameter (named lambda in the article) can be a constant to estimate
//
// Note that alpha > 2 is required to have the existence of a variance!

data
{
	// Size integer
	int<lower = 1> N; // Size of response var

	// Vector data
	vector[N] Y; // Response var
	vector[N] dbh; // Explanatory variable
}

parameters // IMPORTANT: it worth adding constraints, at least to respect the priors, otherwise, a lot of divergence!
{
	real<lower = 0.00000001> alpha_0; // 'intercept' for alpha shifted by 2 (see comments and model)
	real<lower = 0.00000001> alpha_1; // 'slope' for alpha (see comments)
	// real<lower = 0.00000001> lambda_0;
	// real<lower = 0.00000001> lambda_1;
}

model
{
	// Declare variables
	vector[N] alpha;
	// vector[N] lambda;
	for (i in 1:N)
	{
		alpha[i] = alpha_0*dbh[i]^alpha_1; // alpha MUST be an increasing function of dbh, to decrease the average and variance
		// lambda[i] = lambda_0/dbh[i]^lambda_1; // lambda I have no idea
		// alpha[i] = exp(-alpha_1*dbh[i]); // alpha MUST be an increasing function of dbh, to decrease the average
		// lambda[i] = exp(-(dbh[i] - lambda_0)^2/lambda_1^2);
		// lambda[i] = lambda_0*exp(-exp(-lambda_1*dbh[i])); // lambda MUST be an increasing function of dbh
		// alpha[i] = alpha_1*dbh[i]; // alpha MUST be an increasing function of dbh
		// alpha[i] = alpha_0 + alpha_1*dbh[i]; // alpha MUST be an increasing function of dbh
	}

	// Prior
	target += gamma_lpdf(alpha_0 | 1.0/100000, 1.0/100000);
	target += gamma_lpdf(alpha_1 | 1.0/100000, 1.0/100000);
	// target += gamma_lpdf(lambda_0 | 1.0/100000, 1.0/100000);
	// target += gamma_lpdf(lambda_1 | 1.0/100000, 1.0/100000);

	// Model
	for (i in 1:N)
	{
		target += pareto_lpdf(Y[i] | 1, alpha[i]);
		// target += pareto_type_2_lpdf(Y[i] | 1, lambda[i], alpha[i]);
	}
}
