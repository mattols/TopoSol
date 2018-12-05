# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# TOPOSOL VARIABLES
#
# !!! old version of code - will be updated by 01/03/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
require(insol)
require(raster)
require(rgdal)

#create directory
if (dir.exists(paste0(save_path,"TopoSol"))){
  print("Directory exists...")
} else{
  dir.create(paste0(save_path,"TopoSol"))
  print(paste0("Directory created: >", save_path,"TopoSol"))
}
setwd(paste0(save_path,"TopoSol"))

# import DEM
dem0 <- raster(dem_path)

# Read in glacier shapefile
glacier <- readOGR(gsh_path) 

# define months of interest
if (start_month == end_month){
  months = start_month
} else{
  months = c(start_month, end_month)
}
