
#### Aim of prog: Evaluate whether temperature or precipitation of 2050 leads to low growth (< 5th percentile)
## Comments
# To evaluate which environmental variable among temperature and precipitation in 2050 leads to low growth (defined as below the 5th percentile
#	of the growth data for each species), I use factor mapping (Saltelli 2008, chapter 5). More specifically, I use a Monte Carlo Filtering
#	method. It works as follows:
#		1. Define what is behavioural and non-behavioural (here, above or below 5th percentile, respectively)
#		2. Draws the factors from their distributions many times (here, temperature and precipitations distribution for 2050)
#		3. Run the model for each draw
#		4. Classify for each output whether it falls in behavioural or non-behavioural
#		5. Map back onto the factors space X, in order to get for each factor X_i the two distributions:
#			5.1. F = [X_i | behavioural result], and
#			5.2. G = [X_i | non-behavioural result].
#		6. Conclude! If the two distributions F and G are similar, then the parameter seems non-influential. However, if they are clearly
#			different, then X_i is important. Here, i = 1 or 2, for temperature or percipitations

