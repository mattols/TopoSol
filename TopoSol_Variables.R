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

# import DEM
dem0 <- raster(dem_path)

# Read in glacier shapefile
glacier <- readOGR(gsh_path) 

# Convert shapefile (if necessary)
if (class(glacier) != "SpatialPolygons" && class(glacier) != "SpatialPolygonsDataFrame"){
  glacier = SpatialPolygons(list(Polygons(list(Polygon(coordinates(glacier))), "glacier")),proj4string=crs(dem0))
  print("...converted shapefile to SpatialPolygons object")
}

# define months of interest
if (start_month == end_month){
  months = start_month
} else{
  months = c(start_month, end_month)
}
