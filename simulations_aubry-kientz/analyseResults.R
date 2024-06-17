
#### Aim of prog: Plot the posterior of the parameters given delta t

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(tikzDevice)
library(cmdstanr)
library(stringi)

#### Copmute quantiles
## Common variables
indiv = 400
plot = 1

# ls_files = data.table(filename = list.files(path = "./", pattern = paste0("^sigmaProc=0.92_indiv=", indiv, "_plot=", plot, ".*.rds")))
ls_files = data.table(filename = list.files(path = "./", pattern = paste0("^indiv=", indiv, "_plot=", plot, ".*.rds")))
ls_files[, delta_t := as.integer(stri_sub(str = filename, from = stri_locate(str = filename, regex = "deltaT=")[, "end"] + 1,
	to = stri_locate(str = filename, regex = "_results")[, "start"] - 1))]
setorder(ls_files, delta_t)

values = c(beta0 = -4, beta1 = 0.42, beta2 = -0.06, beta3 = 0.076, beta4 = -0.01, sigmaProc = 0.42)
latex_labs = c(beta0 = "\\( \\beta_0 \\)", beta1 = "\\( \\beta_1 \\)", beta2 = "\\( \\beta_2 \\)", beta3 = "\\( \\beta_3 \\)",
	beta4 = "\\( \\beta_4 \\)", sigmaProc = "\\( sigmaproc \\)") #! Change manually sigmaproc to \sigmaproc in the output
ls_parameters = names(values)

## Load all the results and compute quantiles (optimised as the bottleneck is loading the results, so delta_t loop first!)
ls_quantiles = vector(mode = "list", length = ls_files[, .N])
names(ls_quantiles) = ls_files[, delta_t]
for (i in seq_len(ls_files[, .N]))
{
	results = readRDS(ls_files[i, filename])
	quantiles_dt = data.table(parameter = ls_parameters, min = 0, q05 = 0, med = 0, q95 = 0, max = 0, key = "parameter")
	for (param in ls_parameters)
	{
		draws = results$draws(param)
		quantiles_dt[param, c("min", "q05", "med", "q95", "max") := as.list(quantile(x = draws, c(0, 0.025, 0.5, 0.975, 1)))]
	}
	ls_quantiles[[as.character(ls_files[i, delta_t])]] = quantiles_dt
}

ls_quantiles = rbindlist(l = ls_quantiles, idcol = "delta_t")
ls_quantiles[, delta_t := as.numeric(delta_t)]
setkey(ls_quantiles, parameter, delta_t)

#### Plot
for (param in ls_parameters)
{
	tikz(filename = paste0("./", param, "_indiv=", indiv, "_plot=", plot, ".tex"), height = 3, width = 3)
	par(mar = c(5, 5.1, 0.2, 0.2) + 0.1)
	plot(0, type = "n", ann = FALSE, axes = FALSE, xlim = c(1, ls_files[, .N]),
		ylim = c(min(ls_quantiles[param, min(q05)], values[param]), max(ls_quantiles[param, max(q95)], values[param])))
	for (i in seq_len(ls_files[, .N]))
	{
		segments(x0 = i, y0 = ls_quantiles[.(param, ls_files[i, delta_t]), q05], y1 = ls_quantiles[.(param, ls_files[i, delta_t]), q95],
			lwd = 1.25, col = ifelse(ls_files[i, delta_t] > 3, "#CD1A21", "#000000"))
		segments(x0 = i - 0.05, y0 = ls_quantiles[.(param, ls_files[i, delta_t]), q05], x1 = i + 0.05,
			col = ifelse(ls_files[i, delta_t] > 3, "#CD1A21", "#000000"))
		segments(x0 = i - 0.05, y0 = ls_quantiles[.(param, ls_files[i, delta_t]), q95], x1 = i + 0.05,
			col = ifelse(ls_files[i, delta_t] > 3, "#CD1A21", "#000000"))
		points(x = i, y = ls_quantiles[.(param, ls_files[i, delta_t]), med], col = "#295384", bg = "#295384", pch = 23)
	}
	axis(side = 1, at = seq_len(ls_files[, .N]), labels = ls_files[, delta_t])
	title(xlab = "\\( \\Delta t \\)")
	axis(side = 2, las = 1)
	title(ylab = latex_labs[param], line = 4)
	abline(h = values[param], lty = "dashed", col = "#919191", lwd = 0.5)
	dev.off()
}
