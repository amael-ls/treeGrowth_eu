# treeGrowth_eu

## Important notes for diagnosis

Source: https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup

1. **Specifically for low BFMI warnings:**
	- Look at the `pairs` plot to see which primitive parameters are correlated with the `energy__` margin.
	- The primitive parameters that are correlated with the `energy__` margin in the `pairs` plot are a good place to start thinking about reparameterizations. There should be be a negative relationship between `lp__` and `energy__` in the `pairs` plot, but this is not a concern because `lp__` is the logarithm of the posterior kernel rather than a primitive parameter.
2. **Specifically for Rhat, ESS, low BFMI warnings:** You might try setting a higher number of warmup or sampling iterations. Increasing the number of iterations is rarely helpful for resolving divergences/max treedepth warnings.
	- Look at change in bulk-ESS and tail-ESS when the number of iterations increase. If R-hat is less than 1.01 and ESS grows linearly with the number of iterations and eventually exceeds the recommended limit, the mixing is sufficient but MCMC has high autocorrelation requiring a large number of iterations.



## Run order for scripts

1. `prepare_growth_data.R`