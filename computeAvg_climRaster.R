
#### Aim of prog: Compute the climate average of the species distrib climate rasters

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(terra)

#### Compute and save raster
species = as.character(commandArgs(trailingOnly = TRUE))

if (species == "Betula pendula")
	terraOptions(verbose = TRUE, steps = 10)

print(paste0("Running for species <", species, ">"))

ls_var = c("tas", "pr")

for (currentVariable in ls_var)
{
	rs_filename = paste0("./", species, "/", currentVariable, "_mean.tif")
	if (file.exists(rs_filename))
	{
		rs = rast(rs_filename)
		print(paste0("The file <", rs_filename, "> already exists and was loaded"))
	} else {
		rs = rast(paste0("./", species, "/", currentVariable, ".tif"))
		rs_mean = mean(rs)
		writeRaster(x = rs_mean, filename = rs_filename, filetype = "GTiff")
		print(paste("Variable", currentVariable, "done."))
	}
	qt = quantile(values(rs), na.rm = TRUE, probs = c(0, 0.01, 0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975, 0.99, 1))
	saveRDS(qt, paste0("./", species, "/", currentVariable, "_qt.rds"))
}
