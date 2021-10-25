

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
	vector[N] dbh0; // diameter at time t (explanatory variable)
	vector[N] dbh1; // diameter at time t + 1 (response var)
}

parameters // IMPORTANT: it worth adding constraints, at least to respect the priors, otherwise, a lot of divergence!
{
	real alpha;
	real beta;
	real<lower = 0.000001> sigma_eps;
}

model
{
	// Prior
	target += normal_lpdf(alpha | 0, 10000);
	target += normal_lpdf(beta | 0, 10000);
	target += gamma_lpdf(sigma_eps | 1.0/100000, 1.0/100000);

	// Model
	for (i in 1:N)
	{
		target += normal_lpdf(dbh1[i] | alpha + beta*dbh0[i], sigma_eps);
	}
}
