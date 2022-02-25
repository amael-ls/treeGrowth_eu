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

See +@fig:residuals_parent_test1 for the residuals of the `parents' (pa), and +@fig:residuals_child_test1 for the distribution of the residuals of the children. As we can see, the latent states of both parents and children match with the observations. This conclude the first test.

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

## To do list
What I should do (list not necessarily ordered):

- [x] Check the indices
- [x] create the DHARMa data using stan rather than R
- [ ] Check how the residuals are computed (maybe an error there?)
- [x] Plot the likelihood surface (heat map) in the plan (potentialGrowth x dbh_slope)
- [ ] Plot few trees and there states all in the same plot (one tree after each other, create fake years for that)
- [ ] https://esdac.jrc.ec.europa.eu/content/european-soil-database-v20-vector-and-attribute-data
- [ ] Change prior target += gamma_lpdf(latent_dbh[parentsObs_index] | 1.582318^2/10.0, 1.582318/10.0); Make it wring too small. Second option: Make it uniform 0, 2000
- [x] Remove process error: does the latent states (from Stan) match with the latent states recorded in R?
- [x] Check this: https://mc-stan.org/cmdstanr/reference/model-method-generate-quantities.html

