# Few notes to understand some problems on estimating growth
When running the script **growth.R** to fit the growth data, the model does not converge (at least the two error terms). Therefore, I create and analyse dummy data using the script **testForGrowth.R** to see what could be the problem. Surprisingly, even if the model does not converge for the error terms, it is able to find the other parameters, including the latent states. That being said, two parameters (potential growth and dbh slope) are biased. Moreover, I systematically overestimate the first measurement, and underestimate the second measurement for each tree. Obviously the residuals (computed with DHARMa) which compute the differences between the dbh data and the corresponding latent states are not good. To investigate I ran few tests described below.

## Test 1: Fix both error terms
1. I fix both error terms to their real values.
2. I got biased values for potential growth and dbh slope. The slopes for precipitation and the latent states are correctly predicted.
3. I saved the results under the name 'test_errorFixed.rds'

## Test 2: Same as 1 + provide the latent states
1. Both error terms are fixed, and the 'latent' states are known (so it is not latent any more by definition). The likelihood is therefore:
$$ [\theta_G | \varphi, \sigma_p] \propto \prod_i \prod_j [\varphi_{j + 1}^i | \varphi_j^i, \sigma_p, \theta_G] \times [\theta_G] $$
where $\theta_G$ is the vector of parameters for growth, $\sigma_p$ the process error (fixed), and $\varphi$ is the vector of latent states. TO FINISH

## Test 3: Same as 1 + ???

## To do list
What I should do (list not necessarily ordered):
1. Check the indices. Maybe use some prints, create the DHARMa data using stan rather than R, shuffle the data on R
2. Check how the residuals are computed (maybe an error there?)
3. Set one parameter to its true value ()
4. Plot the likelihood surface (heat map) in the plan (potentialGrowth x dbh_slope) which are highly correlated
5. Cannot remember, need to check at the office
