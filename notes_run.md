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

## Abies grandis, file growth-run=1-2022-05-02_21h28.rds

### Print results

```C++
         variable      mean   median     sd    mad       q5      q95 rhat ess_bulk ess_tail
 averageGrowth_mu     -3.37    -3.37   0.07   0.07    -3.48    -3.26 1.09       29      145
 averageGrowth_sd      0.38     0.38   0.02   0.02     0.34     0.42 1.00      409      971
 dbh_slope             0.33     0.33   0.02   0.02     0.31     0.36 1.01      177      361
 pr_slope              0.00     0.00   0.05   0.04    -0.08     0.07 1.08       33      201
 pr_slope2            -0.10    -0.11   0.04   0.04    -0.16    -0.03 1.44        6       28
 tas_slope            -0.19    -0.19   0.03   0.03    -0.24    -0.14 1.02      176      540
 tas_slope2           -0.03    -0.03   0.02   0.02    -0.06    -0.01 1.02      160      487
 ph_slope             -0.08    -0.08   0.04   0.04    -0.14    -0.02 1.04      126      590
 ph_slope2            -0.04    -0.04   0.04   0.04    -0.10     0.02 1.00      286      998
 competition_slope    -0.24    -0.24   0.03   0.03    -0.28    -0.19 1.01      229      826
 sigmaObs[1]           0.02     0.02   0.00   0.00     0.01     0.02 1.05       71      165
 sigmaObs[2]           0.01     0.01   0.00   0.00     0.01     0.02 1.00      383      629
 etaObs[1]             0.15     0.15   0.02   0.02     0.12     0.18 1.00     4556     1940
 etaObs[2]             0.17     0.17   0.02   0.02     0.14     0.20 1.00     4891     2420
 proba[1]              0.01     0.01   0.00   0.00     0.01     0.01 1.00     2523     2320
 proba[2]              0.02     0.02   0.00   0.00     0.02     0.02 1.00     4058     2051
 sigmaProc             0.00     0.00   0.00   0.00     0.00     0.00 1.03      139      524
```

### Messages from pair plots (correlations with energy)

```C++
Correlation for pr_slope2:   -0.202
Correlation for sigmaProc:    0.41
Correlation for sigmaObs[1]:  0.524
```

### Messages about routine and extreme errors

```C++
"France:"
The routine error is 2.227 mm.
The extreme error is 21.204 mm. It occurs with a probability p = 0.0097
Correlation routine <---> extreme =  -0.036

"Germany:"
The routine error is 1.886 mm.
The extreme error is 23.96 mm. It occurs with a probability p = 0.0201
Correlation routine <---> extreme =  -0.001
```

### Posteriors

![Posterior extreme error](./Abies%20grandis/sigmaObs_posterior_1.pdf "Routine error")
![Posterior extreme error](./Abies%20grandis/etaObs_posterior_1.pdf "Extreme error")
![Posterior extreme error](./Abies%20grandis/proba_posterior_1.pdf "Probability")

## Abies grandis, file growth-run=1-2022-05-13_05h41.rds

This run was done with the growth.stan from 2022-05-11, **i.e., with the new structure of random effect, the proba prior is narrowed around 1% ± 0.5**.

### Print results

```C++
 lp__              -6379.58 -6380.16 209.61 206.53 -6721.00 -6040.22 1.03      125      283
 averageGrowth_mu     -3.40    -3.40   0.06   0.06    -3.50    -3.30 1.03      186      420
 averageGrowth_sd      0.38     0.38   0.02   0.02     0.34     0.42 1.01      365      726
 competition_slope    -0.24    -0.24   0.03   0.03    -0.29    -0.19 1.05       69      392
 dbh_slope             0.33     0.33   0.02   0.02     0.31     0.36 1.00      376      799
 ph_slope             -0.09    -0.09   0.04   0.04    -0.15    -0.03 1.01      279      646
 ph_slope2            -0.04    -0.04   0.04   0.04    -0.10     0.02 1.01      375      703
 pr_slope             -0.03    -0.03   0.04   0.04    -0.09     0.04 1.01      146      383
 pr_slope2            -0.08    -0.08   0.03   0.04    -0.13    -0.01 1.10       40       93
 tas_slope            -0.19    -0.19   0.03   0.03    -0.24    -0.14 1.02      173      400
 tas_slope2           -0.03    -0.03   0.02   0.02    -0.06     0.00 1.01      216      506
 sigmaObs[1]           0.01     0.01   0.00   0.00     0.01     0.02 1.05       96      284
 sigmaObs[2]           0.01     0.01   0.00   0.00     0.01     0.01 1.01      316      686
 etaObs[1]             0.18     0.18   0.02   0.02     0.15     0.20 1.00     3455     2446
 etaObs[2]             0.18     0.18   0.02   0.02     0.15     0.21 1.00     3575     2341
 proba[1]              0.00     0.00   0.00   0.00     0.00     0.00 1.00     4181     2126
 proba[2]              0.00     0.00   0.00   0.00     0.00     0.01 1.00     2978     1752
 sigmaProc             0.00     0.00   0.00   0.00     0.00     0.00 1.01      261      568
```





## Abies grandis, file growth-run=1-2022-05-18_00h19.rds

This run was done with the growth.stan from 2022-05-17, **after correction of bug 7 (i.e. using `log_mix` function). We also kept the centred random effect structure but the proba prior is back to Rüger et al 2011**.

### Print results

```c++
 lp__              -5845.62 -5852.23 230.35 237.28 -6211.86 -5446.28 1.06       70      157
 averageGrowth_mu     -3.42    -3.42   0.07   0.06    -3.52    -3.31 1.02      167      437
 averageGrowth_sd      0.38     0.38   0.02   0.03     0.34     0.42 1.06       57      743
 competition_slope    -0.23    -0.23   0.03   0.03    -0.28    -0.19 1.01      273      662
 dbh_slope             0.33     0.33   0.02   0.02     0.30     0.36 1.03      200      608
 ph_slope             -0.09    -0.09   0.04   0.04    -0.15    -0.03 1.01      320      761
 ph_slope2            -0.05    -0.05   0.04   0.04    -0.11     0.02 1.01      295      917
 pr_slope             -0.03    -0.04   0.04   0.04    -0.10     0.04 1.04       50      195
 pr_slope2            -0.05    -0.05   0.04   0.04    -0.11     0.01 1.10       25      108
 tas_slope            -0.19    -0.19   0.03   0.03    -0.25    -0.14 1.04       90      564
 tas_slope2           -0.04    -0.04   0.02   0.02    -0.07     0.00 1.01      255      630
 sigmaObs[1]           0.01     0.01   0.00   0.00     0.01     0.01 1.06       45      100
 sigmaObs[2]           0.01     0.01   0.00   0.00     0.01     0.01 1.01      289      623
 etaObs[1]             0.20     0.20   0.02   0.02     0.17     0.23 1.02      247     1463
 etaObs[2]             0.19     0.19   0.02   0.02     0.16     0.22 1.00     1006     2320
 proba[1]              0.03     0.03   0.00   0.00     0.03     0.04 1.04       65      685
 proba[2]              0.03     0.03   0.00   0.00     0.02     0.04 1.00      878     1454
 sigmaProc             0.00     0.00   0.00   0.00     0.00     0.00 1.03       85      178
```



```C++
Checking sampler transitions treedepth.
Treedepth satisfactory for all transitions.

Checking sampler transitions for divergences.
2 of 7500 (0.03%) transitions ended with a divergence.
These divergent transitions indicate that HMC is not fully able to explore the posterior distribution.
Try increasing adapt delta closer to 1.
If this doesn't remove all divergences, try to reparameterize the model.

Checking E-BFMI - sampler transitions HMC potential energy.
The E-BFMI, 0.15, is below the nominal threshold of 0.30 which suggests that HMC may have trouble exploring the target distribution.
If possible, try to reparameterize the model.

The following parameters had fewer than 0.001 effective draws per transition:
  latent_dbh_parents[15], latent_dbh_parents[159], latent_dbh_parents[366], latent_dbh_parents[388], latent_dbh_parents[410], latent_dbh_parents[653]
Such low values indicate that the effective sample size estimators may be biased high and actual performance may be substantially lower than quoted.

The following parameters had split R-hat greater than 1.05:
  plotEffect[38], plotEffect[39], plotEffect[199], plotEffect[247], pr_slope2, averageGrowth_sd, sigmaObs[1], latent_dbh_parents[4], latent_dbh_parents[5], latent_dbh_parents[15], latent_dbh_parents[23], latent_dbh_parents[43], latent_dbh_parents[60], latent_dbh_parents[116], latent_dbh_parents[117], latent_dbh_parents[127], latent_dbh_parents[159], latent_dbh_parents[161], latent_dbh_parents[167], latent_dbh_parents[209], latent_dbh_parents[221], latent_dbh_parents[241], latent_dbh_parents[252], latent_dbh_parents[254], latent_dbh_parents[309], latent_dbh_parents[340], latent_dbh_parents[358], latent_dbh_parents[366], latent_dbh_parents[388], latent_dbh_parents[393], latent_dbh_parents[394], latent_dbh_parents[395], latent_dbh_parents[399], latent_dbh_parents[410], latent_dbh_parents[413], latent_dbh_parents[419], latent_dbh_parents[421], latent_dbh_parents[424], latent_dbh_parents[425], latent_dbh_parents[458], latent_dbh_parents[472], latent_dbh_parents[474], latent_dbh_parents[487], latent_dbh_parents[509], latent_dbh_parents[533], latent_dbh_parents[535], latent_dbh_parents[574], latent_dbh_parents[589], latent_dbh_parents[599], latent_dbh_parents[616], latent_dbh_parents[618], latent_dbh_parents[621], latent_dbh_parents[645], latent_dbh_parents[653], latent_dbh_parents[662], latent_dbh_parents[665], latent_dbh_parents[701], latent_dbh_parents[715], latent_dbh_parents[745], latent_dbh_parents[766], latent_dbh_parents[779], latent_dbh_parents[782], latent_dbh_parents[796], latent_dbh_parents[812], latent_dbh_parents[829], latent_dbh_parents[861], latent_dbh_parents[862], latent_dbh_parents[871], latent_dbh_parents[878], latent_dbh_parents[890], latent_dbh_parents[907], latent_dbh_parents[910], latent_dbh_parents[911], latent_dbh_parents[917], latent_dbh_parents[927], latent_dbh_parents[931], latent_dbh_parents[932], latent_dbh_parents[944], latent_dbh_parents[963], latent_dbh_parents[969], latent_dbh_parents[978], latent_dbh_parents[980], latent_dbh_parents[981], latent_dbh_parents[993], latent_dbh_parents[1017], latent_dbh_parents[1024], latent_dbh_parents[1032], latent_dbh_parents[1037], latent_dbh_parents[1044], latent_growth[576], latent_growth[577], latent_growth[578], latent_growth[579], latent_growth[580], latent_growth[581], latent_growth[582], latent_growth[583], latent_growth[584], latent_growth[585], latent_growth[4983], latent_growth[4984], latent_growth[4985], averageGrowth[38], averageGrowth[39], averageGrowth[159], averageGrowth[199], averageGrowth[252]
Such high values indicate incomplete mixing and biased estimation.
You should consider regularizating your model with additional prior information or a more effective parameterization.
```

For the problematic plot, it worth noticing that the plots 38, 39, and 199 have only one individual recorded. The plot 247 has 3 individuals, all of them with fast annual growth (measurement error or change of protocol?) between 15.727 and 27.818 mm/year (**these values are in the 90% percentile!**).



## Acer monspessulanum, file growth-run=1-2022-05-02_17h00.rds

### Print results

```C++
          variable      mean    median     sd    mad        q5       q95 rhat ess_bulk ess_tail
 averageGrowth_mu      -3.89     -3.89   0.06   0.06     -4.00     -3.79 1.04       57      157
 averageGrowth_sd       0.44      0.44   0.02   0.02      0.40      0.47 1.01      128      559
 dbh_slope              0.09      0.09   0.02   0.02      0.07      0.12 1.01      130      334
 pr_slope               0.02      0.01   0.05   0.05     -0.05      0.10 1.07       46       91
 pr_slope2              0.01      0.01   0.02   0.02     -0.02      0.03 1.06       44      104
 tas_slope             -0.01     -0.01   0.03   0.03     -0.07      0.05 1.03      120      294
 tas_slope2            -0.04     -0.04   0.02   0.02     -0.07     -0.01 1.09       46      216
 ph_slope               0.09      0.09   0.07   0.07     -0.02      0.20 1.02      217      616
 ph_slope2             -0.11     -0.11   0.03   0.04     -0.17     -0.05 1.02      190      578
 competition_slope     -0.21     -0.21   0.04   0.05     -0.28     -0.14 1.01      115      347
 sigmaObs[1]            0.05      0.05   0.00   0.00      0.04      0.05 1.03      134      649
 etaObs[1]              0.29      0.29   0.04   0.04      0.23      0.35 1.00     2785     2055
 proba[1]               0.01      0.01   0.00   0.00      0.01      0.01 1.00     2563     2028
 sigmaProc              0.00      0.00   0.00   0.00      0.00      0.00 1.16       13      100
```

### Messages from pair plots (correlations with energy)

```C++
Correlation for dbh_slope:  -0.197
Correlation for tas_slope2:  0.18
Correlation for ph_slope2:  -0.1865
Correlation for sigmaProc:   0.871
Correlation for sigmaObs:   -0.318
```

### Messages about routine and extreme errors

```C++
"France:"
The routine error is 3.055 mm.
The extreme error is 19.344 mm. It occurs with a probability p = 0.0078
Correlation routine <---> extreme =  -0.034
```

### Posteriors

![Posterior extreme error](./Acer%20monspessulanum/sigmaObs_posterior_1.pdf "Routine error")
![Posterior extreme error](./Acer%20monspessulanum/etaObs_posterior_1.pdf "Extreme error")
![Posterior extreme error](./Acer%20monspessulanum/proba_posterior_1.pdf "Probability")


## Acer opalus, file growth-run=1-2022-05-03_08h10.rds

### Print results

```C++
          variable      mean    median     sd    mad        q5       q95 rhat ess_bulk ess_tail
 averageGrowth_mu      -4.04     -4.04   0.06   0.06     -4.14     -3.94 1.10       23       95
 averageGrowth_sd       0.38      0.38   0.02   0.02      0.35      0.41 1.04       96      214
 dbh_slope              0.15      0.15   0.01   0.01      0.13      0.18 1.01      174      352
 pr_slope               0.02      0.02   0.03   0.03     -0.02      0.07 1.01       93      213
 pr_slope2             -0.05     -0.05   0.02   0.02     -0.08     -0.01 1.07       39       76
 tas_slope              0.04      0.04   0.02   0.02      0.01      0.07 1.01      194      537
 tas_slope2            -0.02     -0.02   0.02   0.02     -0.04      0.01 1.01      155      240
 ph_slope              -0.10     -0.11   0.08   0.09     -0.24      0.03 1.05       58      186
 ph_slope2              0.03      0.03   0.04   0.04     -0.04      0.09 1.03       96      628
 competition_slope     -0.15     -0.15   0.02   0.02     -0.19     -0.11 1.01      195      489
 sigmaObs[1]            0.01      0.01   0.00   0.00      0.01      0.01 1.09       38       93
 etaObs[1]              0.24      0.24   0.03   0.03      0.20      0.28 1.00     3951     2172
 proba[1]               0.00      0.00   0.00   0.00      0.00      0.01 1.00     3477     2182
 sigmaProc              0.00      0.00   0.00   0.00      0.00      0.00 1.09       24      149
```

### Message from pair plots (correlations with energy)

```C++
Correlation for sigmaProc: 0.49
Correlation for sigmaObs: 0.405
```

### Messages about routine and extreme erorrs

```C++
"France:"
The routine error is 1.093 mm.
The extreme error is 21.537 mm. It occurs with a probability p = 0.0041
Correlation routine <---> extreme =  -0.036
```

### Posteriors

![Posterior extreme error](./Acer%20opalus/sigmaObs_posterior_1.pdf "Routine error")
![Posterior extreme error](./Acer%20opalus/etaObs_posterior_1.pdf "Extreme error")
![Posterior extreme error](./Acer%20opalus/proba_posterior_1.pdf "Probability")

----------------------------

# Copy paste of print for the test logN versus gamma_lupdf

## Gamma_lupdf (with sigmaProc according to Nadja Rüger 2011---personal comm.): Abies grandis/growth-run=1-2022-06-09_13h18.rds

File: "Abies grandis/growth-run=1-2022-06-09_13h18.rds"

```c++
r$> results$print(c("lp__", "averageGrowth", "dbh_slope", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2",
    ^I"ph_slope", "ph_slope2", "competition_slope", "sigmaObs", "etaObs", "proba", "sigmaProc"), max_rows = 20)
          variable     mean   median     sd    mad       q5      q95 rhat ess_bulk ess_tail
 lp__              -1677.36 -1682.84 353.59 355.73 -2244.02 -1064.24 1.07       36       75
 averageGrowth        -3.25    -3.25   0.03   0.03    -3.30    -3.20 1.04       60      398
 dbh_slope             0.23     0.23   0.01   0.01     0.21     0.25 1.05       86      273
 pr_slope             -0.06    -0.06   0.01   0.01    -0.09    -0.04 1.01      134      355
 pr_slope2            -0.03    -0.02   0.01   0.01    -0.05     0.00 1.14       20       60
 tas_slope            -0.22    -0.22   0.02   0.02    -0.25    -0.19 1.02      120      292
 tas_slope2           -0.09    -0.09   0.02   0.02    -0.12    -0.07 1.06       58      168
 ph_slope             -0.07    -0.07   0.01   0.01    -0.09    -0.05 1.01      151      407
 ph_slope2             0.02     0.02   0.01   0.02     0.00     0.04 1.02      145      486
 competition_slope    -0.28    -0.28   0.01   0.02    -0.30    -0.25 1.02      126      370
 sigmaObs[1]           0.05     0.05   0.00   0.00     0.04     0.05 1.01      246      852
 sigmaObs[2]           0.01     0.01   0.00   0.00     0.01     0.02 1.01      588     1202
 etaObs[1]             0.24     0.24   0.01   0.01     0.23     0.26 1.00     1915     2185
 etaObs[2]             0.29     0.29   0.01   0.01     0.27     0.32 1.01      397     1395
 proba[1]              0.10     0.10   0.01   0.01     0.08     0.11 1.00      684     1879
 proba[2]              0.08     0.08   0.01   0.01     0.07     0.09 1.00     1258     1909
 sigmaProc             0.00     0.00   0.00   0.00     0.00     0.00 1.07       38       80
```

## LogN (with sigmaProc according to Nadja Rüger 2011---personal communication): Abies grandis/growth-run=1-2022-06-09_14h23.rds

File: "Abies grandis/growth-run=1-2022-06-09_14h23.rds"

```C++
r$> results$print(c("lp__", "averageGrowth", "dbh_slope", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2",
    ^I"ph_slope", "ph_slope2", "competition_slope", "sigmaObs", "etaObs", "proba", "sigmaProc"), max_rows = 20)
          variable     mean   median     sd    mad       q5      q95 rhat ess_bulk ess_tail
 lp__              -6762.26 -6767.57 175.32 164.22 -7043.06 -6457.99 1.02      103      176
 averageGrowth        -3.13    -3.13   0.04   0.04    -3.18    -3.07 1.01      635     1264
 dbh_slope             0.22     0.22   0.01   0.01     0.20     0.24 1.02      359     1127
 pr_slope             -0.02    -0.02   0.02   0.02    -0.05     0.01 1.00      333      874
 pr_slope2            -0.03    -0.03   0.02   0.02    -0.06    -0.01 1.01      309      731
 tas_slope            -0.16    -0.16   0.02   0.02    -0.19    -0.13 1.00      562     1222
 tas_slope2           -0.05    -0.05   0.01   0.01    -0.07    -0.03 1.00      588     1324
 ph_slope             -0.04    -0.04   0.02   0.02    -0.07    -0.02 1.01      511     1042
 ph_slope2             0.00     0.00   0.02   0.02    -0.03     0.02 1.00      690     1680
 competition_slope    -0.20    -0.20   0.02   0.02    -0.22    -0.17 1.01      482     1462
 sigmaObs[1]           0.01     0.01   0.00   0.00     0.01     0.02 1.03       71       98
 sigmaObs[2]           0.01     0.01   0.00   0.00     0.01     0.01 1.02      329      642
 etaObs[1]             0.17     0.17   0.02   0.02     0.15     0.20 1.00      656     1773
 etaObs[2]             0.19     0.19   0.02   0.02     0.16     0.22 1.00     1602     1804
 proba[1]              0.03     0.03   0.00   0.00     0.03     0.04 1.01      325      797
 proba[2]              0.03     0.03   0.00   0.00     0.02     0.04 1.02      174     1420
 sigmaProc             0.08     0.08   0.00   0.00     0.07     0.08 1.00      289      827
```



```C++
r$> results$print(c("lp__", "averageGrowth", "dbh_slope", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2",
    ^I"ph_slope", "ph_slope2", "competition_slope", "sigmaObs", "etaObs", "proba", "sigmaProc"), max_rows = 20)
          variable     mean   median     sd    mad       q5      q95 rhat ess_bulk ess_tail
 lp__              -6283.50 -6284.04 170.65 167.91 -6562.27 -6000.76 1.02      140      362
 averageGrowth        -3.21    -3.21   0.03   0.03    -3.26    -3.15 1.00      436     1256
 dbh_slope             0.23     0.23   0.01   0.01     0.21     0.25 1.00      410      930
 pr_slope             -0.04    -0.04   0.02   0.02    -0.07     0.00 1.02      135      283
 pr_slope2            -0.03    -0.03   0.02   0.02    -0.05     0.00 1.01      191      364
 tas_slope            -0.17    -0.17   0.02   0.02    -0.20    -0.14 1.00      586     1198
 tas_slope2           -0.06    -0.06   0.01   0.01    -0.08    -0.04 1.00      466     1082
 ph_slope             -0.05    -0.04   0.01   0.01    -0.07    -0.02 1.01      534     1043
 ph_slope2            -0.01    -0.01   0.01   0.01    -0.03     0.02 1.01      667     1245
 competition_slope    -0.20    -0.20   0.02   0.02    -0.23    -0.18 1.00      448     1296
 sigmaObs[1]           0.02     0.02   0.00   0.00     0.01     0.02 1.03       99      241
 sigmaObs[2]           0.01     0.01   0.00   0.00     0.01     0.01 1.01      298      659
 etaObs[1]             0.18     0.18   0.01   0.01     0.15     0.20 1.00     1201     2066
 etaObs[2]             0.19     0.19   0.02   0.02     0.17     0.22 1.00     1079     2322
 proba[1]              0.04     0.04   0.00   0.00     0.03     0.04 1.01      322     1039
 proba[2]              0.03     0.03   0.00   0.00     0.02     0.04 1.01      575     1652
 sigmaProc             0.06     0.06   0.00   0.00     0.06     0.06 1.00      274     1026
```

# Tilia platyphyllos

## Model with sd_a and sd_b, like in Ruger 2011. I almost did not check that model but worked well

```C++
results = readRDS("Tilia platyphyllos/newSigmaObs_growth-run=1-2022-06-20_05h29_sda_sdb_nadjaRuger2011.rds")
r$> results$print(c("lp__", "averageGrowth", "dbh_slope", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2",
    ^I"ph_slope", "ph_slope2", "competition_slope", "sigmaObs_a", "sigmaObs_b", "etaObs", "proba", "sigmaProc"), max_rows = 20)
          variable      mean    median     sd    mad        q5       q95 rhat ess_bulk ess_tail
 lp__              -17016.43 -17023.00 223.20 217.79 -17379.70 -16634.17 1.01      207      572
 averageGrowth         -4.04     -4.04   0.02   0.02     -4.08     -4.01 1.01      477     1821
 dbh_slope              0.14      0.14   0.01   0.01      0.13      0.15 1.01      718     1632
 pr_slope              -0.04     -0.04   0.01   0.01     -0.06     -0.03 1.00      619     1122
 pr_slope2              0.00      0.00   0.01   0.01     -0.02      0.02 1.01      287      653
 tas_slope              0.00      0.00   0.01   0.01     -0.01      0.02 1.01      384     1084
 tas_slope2             0.01      0.01   0.01   0.01      0.00      0.02 1.01      459      809
 ph_slope              -0.06     -0.06   0.01   0.01     -0.08     -0.03 1.02      302     1148
 ph_slope2              0.00      0.00   0.01   0.01     -0.02      0.02 1.01      601     1586
 competition_slope     -0.15     -0.15   0.01   0.01     -0.17     -0.13 1.00      697     1350
 sigmaObs_a[1]          0.01      0.01   0.00   0.00      0.01      0.01 1.01      342      815
 sigmaObs_b[1]          0.00      0.00   0.00   0.00      0.00      0.00 1.00     2928     1851
 etaObs[1]              0.11      0.11   0.01   0.01      0.09      0.13 1.03       97      249
 proba[1]               0.03      0.03   0.00   0.00      0.02      0.03 1.08       34      291
 sigmaProc              0.04      0.04   0.00   0.00      0.04      0.04 1.02      234     1013
```



## Usual lognormal model

```C++	
r$> results = readRDS("Tilia platyphyllos/growth-run=1-2022-06-19_03h29_newlogLik.rds")
r$> results$print(c("lp__", "averageGrowth", "dbh_slope", "pr_slope", "pr_slope2", "tas_slope", "tas_slope2",
    ^I"ph_slope", "ph_slope2", "competition_slope", "sigmaObs", "etaObs", "proba", "sigmaProc"), max_rows = 20)
          variable      mean    median     sd    mad        q5       q95 rhat ess_bulk ess_tail
 lp__              -14500.40 -14499.75 266.69 256.19 -14941.77 -14064.89 1.03      163      368
 averageGrowth         -4.03     -4.03   0.02   0.02     -4.07     -4.00 1.00      647     1326
 dbh_slope              0.13      0.13   0.01   0.01      0.12      0.14 1.00      793     1471
 pr_slope              -0.04     -0.04   0.01   0.01     -0.06     -0.02 1.00      709     1201
 pr_slope2              0.00      0.00   0.01   0.01     -0.02      0.02 1.01      346      853
 tas_slope              0.01      0.01   0.01   0.01     -0.01      0.03 1.00      714     1491
 tas_slope2             0.01      0.01   0.01   0.01      0.00      0.02 1.00      625     1041
 ph_slope              -0.05     -0.05   0.02   0.02     -0.08     -0.03 1.01      577      888
 ph_slope2              0.00      0.00   0.01   0.01     -0.02      0.02 1.00      581     1254
 competition_slope     -0.15     -0.15   0.01   0.01     -0.17     -0.13 1.00      509     1119
 sigmaObs[1]            0.01      0.01   0.00   0.00      0.01      0.01 1.07       64      220
 etaObs[1]              0.15      0.15   0.02   0.02      0.12      0.18 1.01      124      244
 proba[1]               0.02      0.02   0.00   0.00      0.01      0.02 1.06       38      243
 sigmaProc              0.03      0.03   0.00   0.00      0.03      0.04 1.04       97      465
```

