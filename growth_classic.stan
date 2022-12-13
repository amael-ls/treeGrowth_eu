
/*
	Comments (see Appendix C 'Rescaling' for the scaling):
	Reminder (method of first and second moments, i.e., mean and variance or standard deviation):
		- The gamma distribution of stan uses, in this order, shape (alpha) and rate (beta). It can however be reparametrised
			using mean and var: alpha = mean^2/var, beta = mean/var

		- When the gamma distribution is scaled by a constant C, it is equivalent to multiply the rate by C
		
		- The lognormal distribution of stan uses, in this order, meanlog (mu) and sdlog (sigma). It can however be reparametrised
			using mean and sd: mu = log(mean^2/sqrt(sd^2 + mean^2)), sigma = sqrt(log(sd^2/mean^2 + 1))

		- When the lognoraml distribution is scaled by a constant C, it is equivalent to substract log(C) to the meanlog parameter

	Like C++, BUGS, and R, Stan uses 0 to encode false, and 1 to encode true. https://mc-stan.org/docs/functions-reference/logical-functions.html
*/


