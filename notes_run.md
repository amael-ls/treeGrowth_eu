# Important notes for diagnosis

Source: https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup

1. **Specifically for low BFMI warnings:**
	- Look at the `pairs` plot to see which primitive parameters are correlated with the `energy__` margin.
	- The primitive parameters that are correlated with the `energy__` margin in the `pairs` plot are a good place to start thinking about reparameterizations. There should be be a negative relationship between `lp__` and `energy__` in the `pairs` plot, but this is not a concern because `lp__` is the logarithm of the posterior kernel rather than a primitive parameter.
2. **Specifically for Rhat, ESS, low BFMI warnings:** You might try setting a higher number of warmup or sampling iterations. Increasing the number of iterations is rarely helpful for resolving divergences/max treedepth warnings.
	- Look at change in bulk-ESS and tail-ESS when the number of iterations increase. If R-hat is less than 1.01 and ESS grows linearly with the number of iterations and eventually exceeds the recommended limit, the mixing is sufficient but MCMC has high autocorrelation requiring a large number of iterations.

**The version of the Stan model is from 27th April 2022 (on GitHub)**



Copy-paste few notes in stan program (with the gamma distribution, not the logNormal)

​    *//! REDUCE THE VARIANCE OF sigmaProc ! I THINK THIS IS WHY I CAN HAVE SO MUCH*

​    *//! VARIATION IN THE GROWTH. PUT MORE INFORMATION HERE AND IT IS THE OTHER PARAMS*

​    *//! THAT SHOULD HAVE MORE UNCERTAINTIES. THIS IS ALSO WHY I HAVE SO MANY 'CRAZY'*

​    *//! GROWTH IN THE TRAJECTORIES, I SHOULD NOT ALLOW FOR THAT. REDUCING SIGMAPROC*

​    *//! MIGHT ALSO HELP THE PROBA! I WILL HAVE A BETTER SIGNAL FOR THE CHILDREN*

​    *//! WHEN THEY ARE OBVIOUSLY WRONG (THE POSTERIOR OF PROBA ERROR GIVEN MEASURES)*



# Results from lognormal, Fagus sylvatica

## Classic approach, 8000 individuals

```R
          variable     mean   median     sd    mad       q5      q95 rhat ess_bulk ess_tail
 lp__              52237.99 52239.65 112.49 111.49 52051.79 52418.11 1.02      205      776
 averageGrowth        -4.61    -4.61   0.03   0.03    -4.66    -4.55 1.02      172      851
 dbh_slope             0.52     0.52   0.03   0.03     0.48     0.57 1.04       47      511
 dbh_slope2           -0.08    -0.08   0.01   0.01    -0.09    -0.07 1.05       37      474
 pr_slope             -0.03    -0.03   0.01   0.01    -0.05    -0.02 1.00      779     1617
 pr_slope2            -0.02    -0.02   0.01   0.01    -0.02    -0.01 1.01      431     1359
 tas_slope             0.03     0.03   0.01   0.01     0.02     0.05 1.00     1635     1903
 tas_slope2           -0.02    -0.02   0.00   0.00    -0.02    -0.01 1.00     1208     1464
 ph_slope             -0.02    -0.02   0.01   0.01    -0.03    -0.01 1.00     2140     2149
 ph_slope2            -0.02    -0.02   0.01   0.01    -0.03    -0.01 1.00     2140     2239
 competition_slope    -0.17    -0.17   0.01   0.01    -0.18    -0.15 1.02      357     1945
 etaObs[1]             0.19     0.19   0.02   0.02     0.16     0.23 1.01     1087     1837
 etaObs[2]             0.12     0.12   0.01   0.01     0.11     0.14 1.00      495     1168
 etaObs[3]             0.13     0.13   0.03   0.02     0.09     0.18 1.00     3040     1929
 proba[1]              0.02     0.02   0.00   0.00     0.01     0.02 1.03       70      947
 proba[2]              0.03     0.03   0.00   0.00     0.03     0.04 1.04       66      381
 proba[3]              0.03     0.03   0.00   0.00     0.02     0.04 1.00     2884     2187
 sigmaProc             0.54     0.54   0.01   0.01     0.53     0.55 1.03       80     1511
```



## SSM approach, 8000 individuals

```R
         variable      mean    median     sd    mad           q5       q95 rhat ess_bulk ess_tail
 lp__              -50766.29 -50757.30 584.46 569.84   -51723.71 -49839.90 1.02      170      385
 averageGrowth         -5.24     -5.24   0.03   0.03       -5.29     -5.18 1.01      512     1026
 dbh_slope              0.59      0.59   0.03   0.03        0.54      0.64 1.00      622     1358
 dbh_slope2            -0.09     -0.09   0.01   0.01       -0.10     -0.08 1.00      580     1414
 pr_slope              -0.01     -0.01   0.01   0.01       -0.03      0.02 1.00      838     1565
 pr_slope2             -0.02     -0.02   0.01   0.01       -0.03     -0.01 1.01      596     1280
 tas_slope              0.05      0.05   0.01   0.01        0.04      0.06 1.00     1781     2251
 tas_slope2            -0.01     -0.01   0.01   0.01       -0.02      0.00 1.00     1357     1867
 ph_slope              -0.01     -0.01   0.01   0.01       -0.02      0.00 1.00     2472     2908
 ph_slope2             -0.01     -0.01   0.01   0.01       -0.02      0.00 1.00     2703     2868
 competition_slope     -0.21     -0.21   0.01   0.01       -0.22     -0.19 1.00     1261     2210
 etaObs[1]              0.20      0.20   0.03   0.02        0.16      0.25 1.00     2957     2052
 etaObs[2]              0.11      0.11   0.01   0.01        0.10      0.13 1.01      711     2115
 etaObs[3]              0.13      0.13   0.03   0.03        0.09      0.18 1.00     5536     2083
 proba[1]               0.02      0.02   0.00   0.00        0.01      0.02 1.00     1065     1998
 proba[2]               0.04      0.04   0.00   0.00        0.04      0.05 1.00      298      559
 proba[3]               0.03      0.03   0.00   0.00        0.02      0.04 1.00     4211     1896
 sigmaProc              1.21      1.21   0.01   0.01        1.19      1.23 1.01      165      398
```



## SSM approach, 12 000 individuals

```R
          variable      mean    median     sd    mad         q5       q95 rhat ess_bulk ess_tail
 lp__              -76398.94 -76403.55 750.19 772.51  -77595.70 -75152.74 1.04       73      335
 averageGrowth         -5.22     -5.22   0.03   0.03      -5.26     -5.17 1.02      344     1085
 dbh_slope              0.60      0.60   0.02   0.02       0.56      0.64 1.01      683     1505
 dbh_slope2            -0.09     -0.09   0.01   0.01      -0.10     -0.09 1.01      542     1529
 pr_slope               0.02      0.02   0.01   0.01       0.00      0.04 1.00      987     1429
 pr_slope2             -0.03     -0.03   0.01   0.01      -0.04     -0.02 1.01      656     1171
 tas_slope              0.04      0.04   0.01   0.01       0.03      0.05 1.00     1312     2266
 tas_slope2            -0.02     -0.02   0.00   0.00      -0.02     -0.01 1.00     1196     2109
 ph_slope              -0.02     -0.02   0.01   0.01      -0.03     -0.01 1.00     2026     2521
 ph_slope2             -0.01     -0.01   0.00   0.00      -0.02      0.00 1.00     2095     2255
 competition_slope     -0.20     -0.20   0.01   0.01      -0.21     -0.19 1.00     1324     2313
 etaObs[1]              0.19      0.19   0.02   0.02       0.16      0.23 1.00     2315     2491
 etaObs[2]              0.10      0.10   0.01   0.01       0.09      0.11 1.01      741     1553
 etaObs[3]              0.13      0.12   0.03   0.03       0.09      0.18 1.01      433     1773
 proba[1]               0.01      0.01   0.00   0.00       0.01      0.02 1.00      949     1523
 proba[2]               0.05      0.05   0.01   0.01       0.05      0.06 1.02      156      652
 proba[3]               0.03      0.03   0.00   0.00       0.02      0.03 1.00     3808     2257
 sigmaProc              1.20      1.20   0.01   0.01       1.18      1.21 1.04       69      312
```

