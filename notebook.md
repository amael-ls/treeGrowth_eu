---
title: "Few notes to understand some problems on estimating growth"
author: ALS
geometry: "left = 3cm, right = 3cm, top = 2cm, bottom = 2cm"
output: pdf_document
header-includes:
	- \usepackage{xfrac}
	- \usepackage[format=hang,labelsep=endash,font=small,labelfont=bf]{caption}
	- \usepackage{cleveref}
---

# Description of the problem 
When running the script **growth.R** to fit the growth data, the model does not converge (at least the two error terms). Therefore, I create and analyse dummy data using the script **testForGrowth.R** to see what could be the problem. Surprisingly, even if the model does not converge for the error terms, it is able to find the other parameters, including the latent states, back. That being said, two parameters (potential growth and dbh slope) are biased. Moreover, I systematically overestimate the first measurement, and underestimate the second measurement for each tree. Obviously the residuals (computed with DHARMa) which compute the differences between the dbh data and the corresponding latent states are not good. To investigate I ran few tests described below. I give here the real values of the parameters:

| Variable                | Set 1               | Set 2               |
|-------------------------|---------------------|---------------------|
| Potential growth        | $\log(2.47)$        | $\log(2.47)$        |
| Dbh slope               | $\sfrac{0.11}{135}$ | $\sfrac{0.11}{135}$ |
| Precip slope            | -0.34               | 0                   |
| Precip slope square     | 0.0078              | 0                   |
| Measurement error (std) | $\sqrt{3}$          | $\sqrt{3}$          |
| Process error (var)     | 5.68                | 5.68                |

Table: Sets of parameters. We specify in each test whether a parameter should be estimated or is provided {#tbl:sets}

Note that the parameters are given for a non-standardised dbh. This implies that some parameters returned by Stan models must be transformed. This is the case for instance for `potentialGrowth` (the average growth) and `dbh_slope`:
$$
\begin{aligned}
	\text{potentialGrowth}_s &= \text{potentialGrowth}_r - \log \big( \text{sd}(\text{dbh}) \big) \\
	\text{dbh\_slope}_s &= \text{sd}(\text{dbh}) \times \text{dbh\_slope}_r,
\end{aligned}
$$ {#eq:transform}
where the indices r and s are for R language and Stan language respectively. The values from R correspond to the values defined in +@tbl:sets.

The complete model is:
$$
\begin{aligned}
	[\varphi, \theta_G, & \sigma_{\text{obs}}, \sigma_{\text{proc}} | \phi] \propto [\phi_{\text{pa}} | \varphi_{\text{pa}}] \cdot [\phi_{\text{ch}} | \varphi_{\text{ch}}] \times{} \\
	& \prod_{i}^{\text{indiv}}\prod_{j}^{\text{yr}(i)} [\varphi_{j + 1}^{i} | \varphi_{j}^{i}, \theta_G, \sigma_{\text{proc}}] \times{} \\
	& [\varphi_{\text{pa}}] \cdot [\sigma_{\text{obs}}] \cdot [\sigma_{\text{proc}}],
\end{aligned}
$$ {#eq:fullModel}
and will be simplified during the following tests. In equation (@eq:fullModel), the first line is the likelihood, in which the first measurement `parent' (pa) of each tree is compared to the corresponding latent state, and similarly with the following measurements (children, ch). The second line is the Markov process, relating a state to its previous state and the expected growth, and the third line contains the priors. Note that $\varphi_{\text{pa}}$ corresponds to the diffuse initialisation of the first measurement states.

# Test 1: Simple model with no process error
In this test, I removed the process error, and fixed the measurement error and the precipitation slopes (see file `growth_noProcError.stan`). Therefore, the likelihood is defined as follow:
$$
	[\theta_G, \varphi | \phi, \sigma_{\text{obs}}] \propto [\phi_{\text{pa}} | \varphi_{\text{pa}}, \sigma_{\text{obs}}] \cdot [\phi_{\text{ch}} | \varphi_{\text{ch}}, \sigma_{\text{obs}}] \times
		[\varphi_{\text{pa}}] \cdot [\theta_G],
$$ {#eq:modelTest1}
where $\theta_G$ contains `potentialGrowth` (the average growth) and `dbh_slope`, which are the only two parameters to estimate in addition to the states. See set 2 in +@tbl:sets for the parameters. See +@fig:loglik_test1 for the log likelihood.

![Log-likelihood for test 1. Note that I had to transform the estimated parameters provided by Stan, using +@eq:transform](./Tilia_platyphyllos/test_noProcError_measureErrorFixed_prSlopesNull.pdf "Log-likelihood for test 1"){#fig:loglik_test1}

See +@fig:residuals_parent_test1 for the residuals of the `parents' (pa), and +@fig:residuals_child_test1 for the distribution of the residuals of the children.

![Residuals of the parent, using model @eq:modelTest1](./residuals_parent_test1.pdf "Residuals for test 1 (parents)"){#fig:residuals_parent_test1}

![Residuals of the children, using model @eq:modelTest1](./residuals_child_test1.pdf "Residuals for test 1 (children)"){#fig:residuals_child_test1}


```r
sd(residuals) = 1.719633
sqrt(3) = 1.732051 # The observation error
```

# Test 2: as test 1 + climate
This test is done to check that there is no mistake with climate indexing

# Test 3: as test 1 + observation error to estimate
This test is done to check that I can estimate the observation error correctly, in absence of process error

# Test 4: as test 3 + process error given
This test is done to check that I can estimate the observation error correctly, in presence of process error (provided)

<!-- ## Test 1: Fix both error terms
1. I fix both error terms in Stan to their corresponding values from R (see +@tbl:set1).
2. I got biased values for potential growth and dbh slope. The slopes for precipitation are found back, and the latent states have reasonable values although biased.
3. I saved the results under the name 'test_errorFixed.rds'

Added remark after I made test 2: The seed used in this test is the woodstock seed `1969-08-18` (which incidentally equals 1943).

## Test 2: Same as 1 + provide the latent states
Both error terms are fixed, and the 'latent' states are known (so it is not latent any more by definition). The likelihood is therefore:
$$
[\theta_G | \varphi, \sigma_p] \propto \prod_i \prod_j [\varphi_{j + 1}^i | \varphi_j^i, \sigma_p, \theta_G] \times [\theta_G]
$$
where $\theta_G$ is the vector of parameters for growth, $\sigma_p$ the process error (fixed), and $\varphi$ is the vector of latent states (known). The complex index game were checked for the expected_latent_dbh: `count + j - 1` and `climate_index[i] + j - 2`.

## Test 3: Same as 2, with only two parameters to estimate
To test for regression dilution, I fixed the precipitation slopes to 0 and fit the data (generated with `set.seed(123)`). In this case, it worked, I could find back the two remaining parameters to estimate, which are `potentialGrowth` and `dbh_slope`. The likelihood is the same as in test 2, except that $\theta_G$ contains only two parameters.

![Log-likelihood for test 3](./Tilia_platyphyllos/test_errorsFixed_latentGiven_prSlopesNull_seed=123.pdf "Log-likelihood for test 3")

## Test 4: Same as 3, but the latent states are not provided
So this test is like the third one, but the states are not given and must be estimated. Although I am able to estimate the parameters back, there is still a problem in the residuals. Here are the estimated parameters (Tab. @tbl:values_test3), which can be compared to the real values in Tab. @tbl:set1. 

| Variable		   | mean      | sd
|------------------|-----------|---
| lp__  	       |  56446.97 | 90.73
| potential growth | -2.26     | 0.00
| dbh slope  	   |  0.11 	   | 0.00

Table: Values of the estimated parameters for test 3. {#tbl:values_test3}

Plotting the difference in Fig  between measured and simulated (from the chains) shows that the first measurement (parent) and the last measurement (child) should be separated in the analysis of the residuals. Here, each point represent the average residuals over 3000 draws (from Stan). As can be seen, many trees are biased! Ideally, the average would be 0 for each tree.

![Average difference of sampled states versus measured states for the parents (left to the red line) and children (right to the red line)](./meanDiff.pdf "Average difference of simulated versus measured"){#fig:avg_residuals3 width=65%}
 -->

## To do list
What I should do (list not necessarily ordered):

- [x] Check the indices
- [ ] create the DHARMa data using stan rather than R
- [ ] Check how the residuals are computed (maybe an error there?)
- [x] Plot the likelihood surface (heat map) in the plan (potentialGrowth x dbh_slope)
- [ ] Plot few trees and there states all in the same plot (one tree after each other, create fake years for that)
- [ ] https://esdac.jrc.ec.europa.eu/content/european-soil-database-v20-vector-and-attribute-data
- [ ] Change prior target += gamma_lpdf(latent_dbh[parentsObs_index] | 1.582318^2/10.0, 1.582318/10.0); Make it wring too small. Second option: Make it uniform 0, 2000
- [ ] Remove process error: does the latent states (from Stan) match with the latent states recorded in R?
- [ ] Check this: https://mc-stan.org/cmdstanr/reference/model-method-generate-quantities.html

