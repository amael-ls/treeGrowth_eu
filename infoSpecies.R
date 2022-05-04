
#### Aim of prog: To list some useful informations about each species of the data set

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(stringi)

#### Load data
treeData_folder = "/home/amael/project_ssm/inventories/growth/"
treeData = readRDS(paste0(treeData_folder, "standardised_european_growth_data_reshaped.rds"))