## TopoSol model
### Impacts of topographic shading on direct solar radiation for valley glaciers in complex topography
#### [Link to publication](https://doi.org/10.5194/tc-13-29-2019)

#### Requires
1. **Digital Elevation Model (DEM)**
- should be in lat/lon projection
- input full DEM, larger than shapefile (for cast shadow effect)
2. **Glacier shapefile (.shp)**
- with associated files
3. **Define dates of interest and paths to data**
- located in *Toposol_Main.R*



#### Output
- The result of tf.models (**tf**) is a 4-dimensional Rasterstack.
- Each layer is an individual raster including all pixels located within the glacier shapefile. 

**Values show the average change in direct solar radiation for the given time period due to:**

**tf[[1]]** = Slope and aspect

**tf[[2]]** = Cast shadows (component of topographic shading)

**tf[[3]]** = Shaded relief (component of topographic shading)

**tf[[4]]** = Combined (total impact of all topographic attributes)

*if a save_path is defined, an ESRI compatible geotiff file will be saved within a new folder called TopoSol*

**Note:** Values do not show totals in direct insolation, but a change based on different topographic features.

<br>

### New model! 
#### Impacts of topography, DEM resolution, and different DEM products
#### [Link to new publication](https://doi.org/10.3389/feart.2019.00216)
- Model code available under **TopoRad** repository

*TopoRad* is a model that includes all shortwave componenets: direct, diffuse, and reflected solar radiation. Code is available, but may need to be updated.

