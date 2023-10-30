
#### Aim of prog: Check how much area is covered in our data vs the whole species range
## Comments:
# The chorological maps have been downloaded from:
#	Caudullo, Giovanni and Erik, Welk and Jesús, San-Miguel-Ayanz (2017)
#	“Chorological Maps for the Main European Woody Species.” in Data in Brief 12: 662–66.
#	https://doi.org/10.1016/j.dib.2017.05.007.
#
# I also tried the chorolical maps from:
#	Brus, D. J., Hengeveld, G. M., Walvoort, D. J. J., Goedhart, P. W., Heidema, A. H., Nabuurs, G. J. & Gunia, K. (2012)
#	“Statistical mapping of tree species over Europe” in European Journal of Forest Research, 131(1), 145–157.
#	https://doi.org/10.1007/s10342-011-0513-5
#	but I actually did not used these maps. 
#

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(stringi)
library(terra)

#### Load data
## Common variables
ls_species = c("Betula pendula", "Fagus sylvatica", "Picea abies", "Pinus pinaster", "Pinus sylvestris", "Quercus petraea")
ls_shapefiles = vector(mode = "list", length = length(ls_species))
names(ls_shapefiles) = ls_species

dir_shapefiles = "/Users/mistral/ownCloud/database/chorological_maps/"

world = vect("/Users/mistral/t_ownCloud/database/shapefiles/world/WB_countries_Admin0_10m/WB_countries_Admin0_10m.shp")
ocean = vect("/Users/mistral/t_ownCloud/database/shapefiles/world/ne_50m_ocean/ne_50m_ocean.shp")

if (!all(is.valid(world)))
	world = makeValid(world)

if (!all(is.valid(ocean)))
	ocean = makeValid(ocean)

## Species shapefiles
for (species in ls_species)
{
	species_underscore = stri_replace(str = species, replacement = "_", regex = " ")
	currentDir = paste0(dir_shapefiles, species, "/shapefiles/")
	shapefile_name = list.files(path = currentDir, pattern = paste0("^", species_underscore, "_plg.shp"))

	shp = project(vect(paste0(currentDir, shapefile_name)), world)
	if (!all(is.valid(shp)))
		shp = makeValid(shp)
	
	ls_shapefiles[[species]] = crop(x = shp, y = world)
}

ls_ext = lapply(X = ls_shapefiles, FUN = ext)
min_x = min(sapply(X = ls_ext, FUN = function(sp_ext) {return (min(sp_ext)[1])}))
max_x = max(sapply(X = ls_ext, FUN = function(sp_ext) {return (max(sp_ext)[1])})) + 10
min_y = min(sapply(X = ls_ext, FUN = function(sp_ext) {return (min(sp_ext)[2])})) - 5
max_y = max(sapply(X = ls_ext, FUN = function(sp_ext) {return (max(sp_ext)[2])})) + 10

world = crop(x = world, y = ext(c(min_x, max_x, min_y, max_y), xy = FALSE))
ocean = crop(x = ocean, y = ext(c(min_x, max_x, min_y, max_y), xy = FALSE))

for (species in ls_species)
{
	plot(world, axes = FALSE, main = species)
	plot(ls_shapefiles[[species]], add = TRUE, col = "#CD1A2188")
	plot(ocean, add = TRUE, col = "#0C93A322")
}

#### Other data from EFI
quercus = rast("/Users/mistral/ownCloud/database/EU_TreeMap/QuercusRoburPetraea.tif")

world = project(world, quercus)
ocean = project(ocean, quercus)

plot(quercus, axes = FALSE, legend = FALSE)
plot(world, add = TRUE)

fagus = rast("/Users/mistral/ownCloud/database/EU_TreeMap/FagusSpp.tif")
plot(fagus, axes = FALSE, legend = FALSE)
plot(world, add = TRUE)
