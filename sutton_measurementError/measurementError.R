
#### Aim of prog: Investigate on the measurement error
## Explanations
# Around 950 trees have been measured by two people in Sutton Canada. It was done this way:
#	1. The first person measure many trees (more than a hundred in a row) while the second person takes notes
#	2. Then, the roles are reversed: person 2 measures and person 1 takes notes.
# There should have enough measurements in a row to make sure that the second person is not influenced by the first
# 
## Names column data:
#	- Arbre = tree id number
#	- Esp = species
#	- Etat = state (V = alive, there are only living trees)
#	- Multi = is it a multi-trunk or not. N for No
#	- The last two columns are the values collected by the two measurers. Note that the names of the measurers appear in the data provided
#		by the lab Gravel (Sherbrooke, Canada). To preserve the anonymity, these names were replaced by dbh1_in_mm, dbh2_in_mm
#	- The two persons measured the diameters of trees in mm at the height 1.37m (breast height).

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
# library(cmdstanr)
# library(stringi)

#### Load data
treeData = readRDS("./trees_remeasured.rds")

plot(treeData[, dbh1_in_mm], treeData[, dbh2_in_mm])
