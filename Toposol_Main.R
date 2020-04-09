# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 
# Effects of TOPOGRAPHIC SHADING on Direct Solar Radiation
#   code for https://doi.org/10.5194/tc-13-29-2019
#
# Matt Olson - University of Utah Dept of Geography 10/17/2016
#
# Updated 04/2020 - please let me know if you recieve any errors
#     *Code expects the DEM projection to be lat/lon in order calculate
#     the solar geometry in 'tf.models()' (tested using ASTER GDEM)
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #

###################
## DEFINE INPUTS ##
# DEFINE DIRECTORY PATH (path to folder containing all .R TopoSol files)
dirpath = ""
# DEFINE DEM (full path + .tif) (lat/lon projection)
dem_path = ""
# DEFINE GLACIER SHAPEFILE PATH (include .shp)
gsh_path = ""
# DEFINE SAVE PATH (optional location to save data)
#   if provided, will create new folder "TopoSol", and save output
save_path = NULL
# DEFINE YEAR/MONTHS OF INTEREST 
# (use same number for single month)
year = 2016; start_month = 4; end_month = 9
months = seq(start_month,end_month)
###################


# CALL FUNCTIONS
source(paste0(dirpath,"/TopoSol_Functions.R"))
# GENERATE VARIABLES FOR GLACIER BASIN
source(paste0(dirpath,"/TopoSol_Variables.R"))


# CALCULATE TOPOGRAPHIC FORCING (CREATES RASTER STACK & SAVES)
tf <- tf.models(dem0, glacier, months, save_path)

# Output of tf.models() is a 4-dimensional Rasterstack.
#   Each layer is a raster including pixels located within the glacier
#   shapefile. Values show the average change in direct solar radiation
#   for the given time period due to:
#   tf[[1]] = Slope and aspect
#   tf[[2]] = Cast shadows (component of topographic shading)
#   tf[[3]] = Shaded relief (component of topographic shading)
#   tf[[4]] = Combined (total impact of all topographic attributes)
# *more info in The Cryosphere publication 
