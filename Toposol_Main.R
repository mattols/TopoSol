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
#   if provided, will save output of tf.model to specified folder
save_path = NULL
# DEFINE YEAR/MONTHS OF INTEREST 
# (use same number for single month)
year = 2016; start_month = 4; end_month = 9
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

## PLOT THE RESULTS
# basic plot 1 (combined topography)
plot(tf[[4]], main="Combined change in direct SW due to topography", col=rev(heat.colors(20)), 
     legend.args = list(text = bquote(Delta~W~m^-2), side = 1, font = 2, line = -7.1, adj = 0,cex = 0.9))
plot(glacier,add=T); contour(tf[[5]],add=T); legend("topright","Elevation (m)", lwd=1, cex=0.75, bty='n')
# basic plot 2 (total shading)
plot(tf[[2]]+tf[[3]], main="Change in direct SW due to topographic shading", col=rev(blues9), 
     legend.args = list(text = bquote(Delta~W~m^-2), side = 1, font = 2, line = -7.1, adj = 0,cex = 0.9))
plot(glacier,add=T); contour(tf[[5]],add=T); legend("topright","Elevation (m)", lwd=1, cex=0.75, bty='n')
# basic plot 3 (comparison of SW topographic impacts along elevation)
dft = as.data.frame(tf);dft = dft[complete.cases(dft$DEM),]; cols = c("firebrick", "deepskyblue", "blue","grey10")
plot(smooth.spline(dft[,4]~dft$DEM), ylab = "", main = "Changes in direct SW due to topography",
     xlab = "Glacier elevation (m)", col=cols[4],type='l', lwd = 2)
title(ylab = bquote(Change~'in'~irradiance~Delta~W~m^-2),line=2.2)
for (i in 1:3){lines(smooth.spline(dft[,i]~dft$DEM), col = cols[i], lwd = 2)}
abline(h=0,lty=3); legend("topright",names(dft)[1:4], lwd=2, col=cols, cex=0.75, bty='n')
