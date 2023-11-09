
#### Aim of prog: Compare the covered area by the France/Germany/Sweden data to the whole distribution.
## Comments:
# The chorological maps have been downloaded from:
#	Caudullo, Giovanni and Erik, Welk and Jesús, San-Miguel-Ayanz (2017)
#	“Chorological Maps for the Main European Woody Species.” in Data in Brief 12: 662–66.
#	https://doi.org/10.1016/j.dib.2017.05.007.
#
# The european shapefile is from:
#	https://ec.europa.eu/eurostat/web/gisco/
#
# The world shapefile is from:
#	https://datacatalog.worldbank.org/search/dataset/0038272/World-Bank-Official-Boundaries
#

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(stringi)
library(terra)

#### Compute climate range for each species
## Common variables
ls_species = c("Betula pendula", "Fagus sylvatica", "Picea abies", "Pinus pinaster", "Pinus sylvestris", "Quercus petraea")

dir_shapefiles = "/home/amael/shapefiles/trees/europe/chorological_maps/"

europe = vect("/home/amael/shapefiles/europe/continent/europe.shp")
world = vect("/home/amael/shapefiles/world/WB_countries_Admin0_10m/WB_countries_Admin0_10m.shp")

europe = project(europe, world)

france = europe[europe$NAME_LATN == "France"]
germany = europe[europe$NAME_LATN == "Deutschland"]
sweden = europe[europe$NAME_LATN == "Sverige"]

areaCover = data.table(speciesName_sci = ls_species, totalArea = -Inf, area_fr = -Inf, area_de = -Inf, area_sw = -Inf)
setkey(areaCover, speciesName_sci)

## Load species shapefiles and world climate rasters to compute range
for (species in ls_species)
{
	species_underscore = stri_replace(str = species, replacement = "_", regex = " ")
	currentDir = paste0(dir_shapefiles, species, "/shapefiles/")
	shapefile_name = list.files(path = currentDir, pattern = paste0("^", species_underscore, "_plg.shp"))

	shp = project(vect(paste0(currentDir, shapefile_name)), world)
	if (!all(is.valid(shp)))
		shp = makeValid(shp)
	
	shp = crop(x = shp, y = world)
	
	areaCover[species, totalArea := sum(expanse(x = shp, unit = "km", transform = TRUE))]
	areaCover[species, area_fr := sum(expanse(terra::intersect(x = shp, y = france), unit = "km", transform = TRUE))]
	areaCover[species, area_de := sum(expanse(terra::intersect(x = shp, y = germany), unit = "km", transform = TRUE))]
	areaCover[species, area_sw := sum(expanse(terra::intersect(x = shp, y = sweden), unit = "km", transform = TRUE))]
}

areaCover[, percentCover := (area_fr + area_de + area_sw)*100/totalArea]

saveRDS(areaCover, "./areaCover.rds")
