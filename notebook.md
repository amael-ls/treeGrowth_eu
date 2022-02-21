# Few notes to understand some problems on estimating growth
When running the script **growth.R** to fit the growth data, the model does not converge (at least the two error terms). Therefore, I create and analyse dummy data using the script **testForGrowth.R** to see what could be the problem. Surprisingly, even if the model does not converge for the error terms, it is able to find the other parameters, including the latent states. That being said, two parameters (potential growth and dbh slope) are biased. Moreover, I systematically overestimate the first measurement, and underestimate the second measurement for each tree. Obviously the residuals (computed with DHARMa) which compute the differences between the dbh data and the corresponding latent states are not good. To investigate I ran few tests described below.

## Test 1: Fix both error terms
1. I fix both error terms to their real values.
2. I got biased values for potential growth and dbh slope. The slopes for precipitation and the latent states are correctly predicted.
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

![Log-likelihood for test 3. Only two parameters to estimate: `potentialGrowth` and `dbh_slope`](./Tilia_platyphyllos/test_errorsFixed_latentGiven_prSlopesNull_seed=123.pdf "Log-likelihood for test 3")

## Test 4: Same as 3, but the latent states are not provided
So this test is like the third, but the states are not given and must be estimated. So far it works

## To do list
What I should do (list not necessarily ordered):

- [x] Check the indices
- [ ] create the DHARMa data using stan rather than R
- [ ] Check how the residuals are computed (maybe an error there?)
- [x] Plot the likelihood surface (heat map) in the plan (potentialGrowth x dbh_slope)
- [ ] Plot few trees and there states all in the same plot (one tree after each other, create fake years for that)

