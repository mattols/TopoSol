# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 
# Effects of TOPOGRAPHIC SHADING on Direct Solar
#
# Matt Olson - University of Utah Dept of Geography 10/17/2016
#
# !!! old version of code - will be updated by 01/03/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #

###################
## DEFINE INPUTS ##
# DEFINE DIRECTORY PATH (location of all files)
dirpath = ""
# DEFINE DEM (full path + .tif)
dem_path = ""
# DEFINE GLACIER SHAPEFILE PATH (include .shp)
gsh_path = ""
# DEFINE SAVE PATH (creates new folder "TopoSol")
save_path = ""
# DEFINE YEAR/MONTHS OF INTEREST 
# (use same number for single month)
year = 2016; start_month = 4; end_month = 9
###################


# CALL FUNCTIONS
source(paste0(dirpath,"TopoSol_Functions.R"))
# GENERATE VARIABLES FOR GLACIER BASIN
source(paste0(dirpath,"TopoSol_Variables.R"))


# CALCULATE TOPOGRAPHIC FORCING (CREATES RASTER STACK & SAVES)
tf <- tf.models(dem0, glacier, months,)
