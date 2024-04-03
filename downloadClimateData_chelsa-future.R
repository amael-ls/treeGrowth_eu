#### Aim of prog: download Chelsa climate projections for 2050
## Explanations
# 1. This script downloads climate projections for the world.
# 2. The downloaded data are tif.
# 4. Documentation is in English
# 5. To cite the data, check in the doc:
# CHELSA projections should be cited as:
#	Brun, Philipp, Niklaus E. Zimmermann, Chantal Hari, Loïc Pellissier, and Dirk Nikolaus Karger, 2022
#	Global Climate-Related Predictors at Kilometer Resolution for the Past and Future
#	Earth System Science Data 14 (12)
#	https://doi.org/10.5194/essd-14-5573-2022.
# 6. Naming convention for files: CHELSA_<variable>_mon_<model>_<rcp>_[other stuff]_<month>_<time period>.tif
#		where month is a number between 01 and 12, and time period is (for me) 2041-2060
#
# Variables and units:
#	bio1	= Annual Mean Temperature [°C*10]
#	bio12	= Annual Precipitation [mm/year]
#
# Shared socioeconomic pathways:
#	More informations can be found at https://www.dkrz.de/en/communication/climate-simulations/cmip6-en/the-ssp-scenarios
#	SSP1: The sustainable and “green” pathway describes an increasingly sustainable world.
#	SSP3: Regional rivalry. A revival of nationalism and regional conflicts pushes global issues into the background.
#	SSP5: Fossil-fueled Development. Global markets are increasingly integrated, leading to innovations and technological progress.
#	Then, they are combined with RCP scenarios, to get for instance:
#		- SSP585: This scenario represents the upper boundary of the range of scenarios described in the literature
#		- SSP370: This scenario is in the upper-middle part of the full range of scenarios
# 		- SSP126: Remake of the optimistic scenario RCP2.6, designed with the aim of simulating a development that is compatible with
#			the 2°C target. This scenario, too, assumes climate protection measures being taken.
#

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500, timeout = 200)

##################################################################
############            DATA FOR THE WORLD            ############
##################################################################

#### Define common variables
## Folder where data will be downloaded
folderDownloadData = "~/scratch/Chelsa/"

if (!dir.exists(folderDownloadData))
	stop(paste0("Folder '", folderDownloadData, "' does not exist. Create the folder before downloading data"))

## Link to download the data
downloadLink = "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/2041-2070/"

#### Quick and dirty download for monthly timeseries for the whole world
# The pattern link is:
# https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/2041-2070/GFDL-ESM4/ssp370/bio/
#	CHELSA_bio12_2041-2070_gfdl-esm4_ssp370_V.2.1.tif

ls_models = c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll")
ssp = "ssp370" # There is also ssp126, and ssp585

for (var in c("bio1", "bio12"))
{
	for (model in ls_models)
	{
		url_file = paste0(downloadLink, toupper(model), "/", ssp, "/bio/CHELSA_", var, "_2041-2070_", model, "_", ssp, "_V.2.1.tif")
		
		exitDir = paste0(folderDownloadData, var, "/", model, "/")
		if (!dir.exists(exitDir))
			dir.create(path = exitDir, recursive = TRUE)

		saveName = paste0(var, "_", model, ".tif")
		download.file(url = url_file, destfile = paste0(exitDir, saveName))
	}
	print(paste("Variable ", var, "done"))
}
