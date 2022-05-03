# Important notes for diagnosis

Source: https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup

1. **Specifically for low BFMI warnings:**
	- Look at the `pairs` plot to see which primitive parameters are correlated with the `energy__` margin.
	- The primitive parameters that are correlated with the `energy__` margin in the `pairs` plot are a good place to start thinking about reparameterizations. There should be be a negative relationship between `lp__` and `energy__` in the `pairs` plot, but this is not a concern because `lp__` is the logarithm of the posterior kernel rather than a primitive parameter.
2. **Specifically for Rhat, ESS, low BFMI warnings:** You might try setting a higher number of warmup or sampling iterations. Increasing the number of iterations is rarely helpful for resolving divergences/max treedepth warnings.
	- Look at change in bulk-ESS and tail-ESS when the number of iterations increase. If R-hat is less than 1.01 and ESS grows linearly with the number of iterations and eventually exceeds the recommended limit, the mixing is sufficient but MCMC has high autocorrelation requiring a large number of iterations.



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

### Messages from pair plots

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

## Acer monspessulanum, file growth-run=1-2022-05-02_17h00.rds

### Print results



### Messages from pair plots

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

