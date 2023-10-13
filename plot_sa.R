
#### Aim of prog: plot sensitivity analysis, one variable with all the species and methods per panel, for all variables

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(tikzDevice)
library(stringi)

#### Tool functions
## Function to plot sensitivity analysis, one variable per panel with all the species/methods
plot_list_sa = function(sa_ssm_data_list, sa_classic_data_list, ls_data = c("dbh", "pr", "tas", "ph", "ba"),
	ls_titles = c(dbh = "Diameter", pr = "Precipitation", tas = "Temperature", ph = "pH", ba = "Basal area"),
	ext = NULL, relative = TRUE, normalise = TRUE)
{
	# Check data
	if (length(sa_ssm_data_list) != length(sa_classic_data_list))
		stop("Lists should be of the same size")

	ls_species = names(sa_ssm_data_list)
	if (!all(ls_species %in% names(sa_classic_data_list)))
		stop("Lists should contain the same species names")

	if (length(ls_titles) != length(ls_data))
		stop("ls_data and ls_titles should be of the same size")

	if (!all(ls_data %in% sapply(sa_ssm_data_list, function(zz) {return(unique(zz$sa[, parameters]))})))
		stop("Some data are not in the SSM list")
	
	if (!all(ls_data %in% sapply(sa_classic_data_list, function(zz) {return(unique(zz$sa[, parameters]))})))
		stop("Some data are not in the classic list")

	if (!is.null(ext))
	{
		if (length(ext) > 1)
		{
			ext = ext[1]
			warning("Only first extension was kept")
		}
		if (!(ext %in% c(".pdf", ".tex")))
			stop("File type not recognised. Only .pdf, and .tex are recognised")
	}

	# Common variables
	ls_species_label = stri_replace_all(str = ls_species, replacement = " ", regex = "_")
	ls_species_label = stri_trans_totitle(ls_species_label, opts_brkiter = stri_opts_brkiter(type = "sentence"))

	# Multiply data by var(y) if not relative
	if (!relative)
	{
		for (currentSpecies in ls_species)
		{
			sa_ssm_data_list[[currentSpecies]]$sa = merge.data.table(sa_ssm_data_list[[currentSpecies]]$sa,
				data.table(run = seq_along(sa_ssm_data_list[[currentSpecies]]$var_y), var_y = sa_ssm_data_list[[currentSpecies]]$var_y),
				by = "run")
			sa_classic_data_list[[currentSpecies]]$sa = merge.data.table(sa_classic_data_list[[currentSpecies]]$sa,
				data.table(run = seq_along(sa_classic_data_list[[currentSpecies]]$var_y), var_y = sa_classic_data_list[[currentSpecies]]$var_y),
				by = "run")

			sa_ssm_data_list[[currentSpecies]]$sa[, original := original*var_y]
			sa_classic_data_list[[currentSpecies]]$sa[, original := original*var_y]
		}
	}

	# Define empty plot
	x_max = 2.25*length(sa_ssm_data_list)
	y_max = max(max(sapply(sa_ssm_data_list, function(zz) {return (max(zz$sa[sensitivity == "Si", original]))})),
		max(sapply(sa_classic_data_list, function(zz) {return (max(zz$sa[sensitivity == "Si", original]))})))
	
	for (currentVar in ls_data)
	{
		if (!is.null(ext))
		{
			if (ext == ".pdf")
				pdf(file = paste0(currentVar, ext), height = 6, width = 9.708204) # Golden ratio

			if (ext == ".tex")
				tikz(file = paste0(currentVar, ext), height = 6, width = 9.708204) # Golden ratio
		}
		plot(0, type = "n", xlim = c(0, x_max), ylim = c(0, y_max), ylab = "SA", main = ls_titles[currentVar],
			xlab = "Species", las = 1, xaxt = "n")

		x_orig = 0.5
		x_pos_label = numeric(length(ls_species))
		names(x_pos_label) = ls_species

		for (currentSpecies in ls_species)
		{
			x_pos_label[currentSpecies] = x_orig
			x1 = x_orig - 0.05
			x2 = x_orig + 0.05
			y1 = quantile(sa_ssm_data_list[[currentSpecies]]$sa[(sensitivity == "Si") & (parameters == currentVar), original],
				c(0.05, 0.5, 0.95))
			y2 = quantile(sa_classic_data_list[[currentSpecies]]$sa[(sensitivity == "Si") & (parameters == currentVar), original],
				c(0.05, 0.5, 0.95))

			segments(x0 = x1, y0 = y1["5%"], x1 = x1, y1 = y1["95%"], col = "#E9851D")
			points(x = x1, y = y1["50%"], pch = 19, cex = 1, col = "#E9851D")
			points(x = x1, y = y1["5%"], pch = "-", cex = 2, col = "#E9851D")
			points(x = x1, y = y1["95%"], pch = "-", cex = 2, col = "#E9851D")

			segments(x0 = x2, y0 = y2["5%"], x1 = x2, y1 = y2["95%"], col = "#2E77AB")
			points(x = x2, y = y2["50%"], pch = 19, cex = 1, col = "#2E77AB")
			points(x = x2, y = y2["5%"], pch = "-", cex = 2, col = "#2E77AB")
			points(x = x2, y = y2["95%"], pch = "-", cex = 2, col = "#2E77AB")

			abline(v = x_orig, lwd = 0.5, lty = "dashed", col = "#55555555")

			x_orig = x_orig + 2.5
		}
		axis(side = 1, at = x_pos_label, labels = ls_species_label)

		legend(x = "topright", legend = c("SSM", "Classic"), lty = 1, lwd = 3, box.lwd = 0, title = "Model", col = c("#E9851D", "#2E77AB"))
		if (!is.null(ext))
			dev.off()
	}
}

## Function to plot proportion sensitivity analysis, one variable per panel with all the species/methods
plot_list_sa_prop = function(sa_ssm_data_list, sa_classic_data_list, ls_data = c("dbh", "pr", "tas", "ph", "ba"),
	legend_text = c("dbh", "precip", "temp", "pH", "B.A."), ext = NULL, relative = TRUE, normalise = TRUE)
{
	# Check data
	if (length(sa_ssm_data_list) != length(sa_classic_data_list))
		stop("Lists should be of the same size")

	ls_species = names(sa_ssm_data_list)
	if (!all(ls_species %in% names(sa_classic_data_list)))
		stop("Lists should contain the same species names")

	n_data = length(ls_data)
	if (n_data > 5)
		warning("The colours might not be acurate as I am using the palette Johnson from MetBrewer")
	
	if (length(legend_text) != n_data)
		stop("ls_data and legend_text should be of the same size")

	if (!all(ls_data %in% sapply(sa_ssm_data_list, function(zz) {return(unique(zz$sa[, parameters]))})))
		stop("Some data are not in the SSM list")
	
	if (!all(ls_data %in% sapply(sa_classic_data_list, function(zz) {return(unique(zz$sa[, parameters]))})))
		stop("Some data are not in the classic list")

	if (!is.null(ext))
	{
		if (length(ext) > 1)
		{
			ext = ext[1]
			warning("Only first extension was kept")
		}
		if (!(ext %in% c(".pdf", ".tex")))
			stop("File type not recognised. Only .pdf, and .tex are recognised")
	}

	# Common variables
	ls_species_label = stri_replace_all(str = ls_species, replacement = " ", regex = "_")
	ls_species_label = stri_trans_totitle(ls_species_label, opts_brkiter = stri_opts_brkiter(type = "sentence"))

	colours = MetBrewer::met.brewer("Johnson", n_data)
	colours_str = grDevices::colorRampPalette(colours)(n_data)
	names(colours_str) = ls_data

	x_max = 2.25*length(sa_ssm_data_list)

	x_orig = 0.5
	x_pos_label = numeric(length(ls_species))
	names(x_pos_label) = ls_species
	
	delta = 0.15

	ls_rect_dt_y1 = vector(mode = "list", length = length(ls_species))
	names(ls_rect_dt_y1) = ls_species
	ls_rect_dt_y2 = vector(mode = "list", length = length(ls_species))
	names(ls_rect_dt_y2) = ls_species

	# Multiply data by var(y) if not relative
	if (!relative)
	{
		for (currentSpecies in ls_species)
		{
			sa_ssm_data_list[[currentSpecies]]$sa = merge.data.table(sa_ssm_data_list[[currentSpecies]]$sa,
				data.table(run = seq_along(sa_ssm_data_list[[currentSpecies]]$var_y), var_y = sa_ssm_data_list[[currentSpecies]]$var_y),
				by = "run")
			sa_classic_data_list[[currentSpecies]]$sa = merge.data.table(sa_classic_data_list[[currentSpecies]]$sa,
				data.table(run = seq_along(sa_classic_data_list[[currentSpecies]]$var_y), var_y = sa_classic_data_list[[currentSpecies]]$var_y),
				by = "run")

			sa_ssm_data_list[[currentSpecies]]$sa[, original := original*var_y]
			sa_classic_data_list[[currentSpecies]]$sa[, original := original*var_y]
		}
	}

	# Compute rectangles
	for (currentSpecies in ls_species)
	{
		level1 = 0
		level2 = 0
		x1 = x_orig - 0.25
		x2 = x_orig + 0.25

		rect_dt_y1 = data.table(variables = ls_data, xleft = x1 - delta, ybottom = numeric(n_data),
			xright = x1 + delta, ytop = numeric(n_data))
		rect_dt_y2 = data.table(variables = ls_data, xleft = x2 - delta, ybottom = numeric(n_data),
			xright = x2 + delta, ytop = numeric(n_data))
		
		setkey(rect_dt_y1, variables)
		setkey(rect_dt_y2, variables)

		for (currentVar in ls_data)
		{
			x_pos_label[currentSpecies] = x_orig

			y1 = quantile(sa_ssm_data_list[[currentSpecies]]$sa[(sensitivity == "Si") & (parameters == currentVar), original],
				c(0.05, 0.5, 0.95))
			y2 = quantile(sa_classic_data_list[[currentSpecies]]$sa[(sensitivity == "Si") & (parameters == currentVar), original],
				c(0.05, 0.5, 0.95))

			# Rectangles
			rect_dt_y1[currentVar, c("ybottom", "ytop") := .(level1, level1 + unname(y1["50%"]))]
			rect_dt_y2[currentVar, c("ybottom", "ytop") := .(level2, level2 + unname(y2["50%"]))]

			level1 = level1 + unname(y1["50%"])
			level2 = level2 + unname(y2["50%"])
		}

		# Normalisation (to get max = 1, the sensobol package have some uncertainties)
		if (normalise)
		{
			rect_dt_y1[, c("ybottom", "ytop") := .(ybottom/level1, ytop/level1)]
			rect_dt_y2[, c("ybottom", "ytop") := .(ybottom/level2, ytop/level2)]
		}

		ls_rect_dt_y1[[currentSpecies]] = rect_dt_y1
		ls_rect_dt_y2[[currentSpecies]] = rect_dt_y2

		x_orig = x_orig + 2.5
	}

	rect_dt_y1 = rbindlist(l = ls_rect_dt_y1, idcol = "species")
	setkey(rect_dt_y1, species, variables)
	rect_dt_y2 = rbindlist(l = ls_rect_dt_y2, idcol = "species")
	setkey(rect_dt_y2, species, variables)

	# Plot rectangles
	if (!is.null(ext))
	{
		if (ext == ".pdf")
			pdf(file = paste0("all_in_one", ext), height = 6, width = 9.708204) # Golden ratio

		if (ext == ".tex")
			tikz(file = paste0("all_in_one", ext), height = 6, width = 9.708204) # Golden ratio
	}

	if (normalise)
	{
		y_max = 1
	} else {
		y_max = max(rect_dt_y1[, sum(ytop), by = .(species, variables)][, max(V1)], rect_dt_y2[, sum(ytop), by = .(species, variables)][, max(V1)])
	}

	plot(0, type = "n", xlim = c(0, x_max), ylim = c(0, y_max), ylab = "SA", main = "", xlab = "Species", las = 1, xaxt = "n")

	for (currentSpecies in ls_species)
	{
		for (currentVar in ls_data)
		{
			rect(xleft = rect_dt_y1[.(currentSpecies, currentVar), xleft], ybottom = rect_dt_y1[.(currentSpecies, currentVar), ybottom],
				xright = rect_dt_y1[.(currentSpecies, currentVar), xright], ytop = rect_dt_y1[.(currentSpecies, currentVar), ytop],
				density = NULL, col = colours_str[currentVar], border = NA)

			rect(xleft = rect_dt_y2[.(currentSpecies, currentVar), xleft], ybottom = rect_dt_y2[.(currentSpecies, currentVar), ybottom],
				xright = rect_dt_y2[.(currentSpecies, currentVar), xright], ytop = rect_dt_y2[.(currentSpecies, currentVar), ytop],
				density = NULL, col = colours_str[currentVar], border = NA)
		}
		# Vertical line species
		abline(v = x_orig, lwd = 0.5, lty = "dashed", col = "#55555555")
	}
	axis(side = 1, at = x_pos_label, labels = ls_species_label)

	legend(x = "topleft", legend = legend_text, box.lwd = 0, fill = colours_str, horiz = TRUE,
		xpd = TRUE, inset = c(0, -0.1))
	if (!is.null(ext))
		dev.off()
}

#### Load data
## Common variables
ls_species = c("Betula pendula", "Fagus sylvatica", "Picea abies", "Pinus pinaster", "Pinus sylvestris", "Quercus petraea")
# ls_species = "Pinus pinaster" #! WATCHOUT, LOADING FOR ONE SPECIES ONLY !!!!

sa_ssm_data_list = vector(mode = "list", length = length(ls_species))
names(sa_ssm_data_list) = ls_species

sa_classic_data_list = vector(mode = "list", length = length(ls_species))
names(sa_classic_data_list) = ls_species

for (species in ls_species)
{
	sa_ssm_data_list[[species]] = readRDS(paste0("./", species, "/sa_ssm_data.rds")) #! WATCHOUT, LOADING FOR n_param = 1 !!!!
	sa_classic_data_list[[species]] = readRDS(paste0("./", species, "/sa_classic_data.rds")) #! WATCHOUT, LOADING FOR n_param = 1 !!!!
}

#### Plot
# plot_list_sa(sa_ssm_data_list = sa_ssm_data_list, sa_classic_data_list = sa_classic_data_list,
# 	ls_data = c("dbh", "precip", "tas", "ph", "standBasalArea"),
# 	ls_titles = c(dbh = "Diameter", precip = "Precipitation", tas = "Temperature", ph = "pH", standBasalArea = "Basal area"),
# 	ext = ".pdf")


# plot_list_sa_prop(sa_ssm_data_list = sa_ssm_data_list, sa_classic_data_list = sa_classic_data_list, ext = ".pdf", relative = TRUE)

pdf("all.pdf", height = 6, width = 9.708204) # Golden ratio
plot_list_sa(sa_ssm_data_list = sa_ssm_data_list, sa_classic_data_list = sa_classic_data_list, relative = TRUE, normalise = TRUE)
plot_list_sa_prop(sa_ssm_data_list = sa_ssm_data_list, sa_classic_data_list = sa_classic_data_list, ext = NULL,
	relative = FALSE, normalise = FALSE)
plot_list_sa_prop(sa_ssm_data_list = sa_ssm_data_list, sa_classic_data_list = sa_classic_data_list, ext = NULL,
	relative = TRUE, normalise = TRUE)
dev.off()

for (currentSpecies in ls_species)
{
	print(paste0("-----   ", currentSpecies, "   -----"))
	print(mean(sa_ssm_data_list[[currentSpecies]]$var_y))
	print(mean(sa_classic_data_list[[currentSpecies]]$var_y))
	# sd_dbh_ssm = readRDS(paste0("./", currentSpecies, "/1_stanData.rds"))$sd_dbh
	# sd_dbh_classic = readRDS(paste0("./", currentSpecies, "/1_stanData_classic.rds"))$sd_dbh
	# print(paste0("-----   ", currentSpecies, "   -----"))
	# print(mean(sd_dbh_ssm*sa_ssm_data_list[[currentSpecies]]$var_y))
	# print(mean(sd_dbh_classic*sa_classic_data_list[[currentSpecies]]$var_y))
}

growth = readRDS(paste0("./", currentSpecies, "/1_stanData.rds"))$avg_yearly_growth_obs
range(growth)
hist(growth)
var(growth)

sa_ssm_data_list = list(sa_ssm_data_list[[2]])
names(sa_ssm_data_list) = ls_species[2]
sa_classic_data_list = list(sa_classic_data_list[[2]])
names(sa_classic_data_list) = ls_species[2]

pdf("sa_relative.pdf", height = 6, width = 9.708204) # Golden ratio
plot_list_sa_prop(sa_ssm_data_list, sa_classic_data_list = sa_classic_data_list, ext = NULL,
	relative = TRUE, normalise = TRUE)
dev.off()

pdf("sa_absolute.pdf", height = 6, width = 9.708204) # Golden ratio
plot_list_sa_prop(sa_ssm_data_list, sa_classic_data_list = sa_classic_data_list, ext = NULL,
	relative = FALSE, normalise = FALSE)
dev.off()
