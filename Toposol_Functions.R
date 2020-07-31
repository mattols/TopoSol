# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 
# TOPOSOL FUNCTIONS
#
# !!! outdated version of code - will be updated by 01/03/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
require(insol)
require(raster)
require(rgdal)

# MAKE RASTER
make.raster <- function(matrix, dem){
  raster(matrix,
         xmn=dem@extent@xmin, xmx=dem@extent@xmax,
         ymn=dem@extent@ymin, ymx=dem@extent@ymax, 
         crs=crs(dem))
}

# CHANGE N-FACING GLACIER FOR "-" ASPECT - ALTERS AZIMUTH ANGLES
change.aspect <- function(AZIMUTH, dir.shift = "S"){
  if (dir.shift == 'W'){
    AZIMUTH = AZIMUTH + 90 # (rotate Aspect 90deg CCW) WEST ASPECT
  } else if (dir.shift == 'S'){
    AZIMUTH = AZIMUTH  + 180 # (rotate Aspect 180deg CCW) SOUTH ASPECT
  } else if (dir.shift == 'E'){
    AZIMUTH = AZIMUTH - 90 #  (rotate Aspect 90deg CW) EAST ASPECT
  }
  
  if (any(AZIMUTH > 360)){
    AZIMUTH[AZIMUTH > 360] = AZIMUTH[AZIMUTH > 360] - 360
  } else if (any(AZIMUTH < 0)){
    AZIMUTH[AZIMUTH < 0] = AZIMUTH[AZIMUTH < 0] + 360
  }
  return(AZIMUTH)
}

# CROP raster based on shape and buffer 
# (KEEP AT 0.5 to incorporate nearby terrain effects)
crop.raster <- function(stk, shp, buffer = 0.05){
  r <- crop(stk,extent(xmin(shp)-buffer,xmax(shp)+buffer,
                       ymin(shp)-buffer,ymax(shp)+buffer))
  return(r)
}

julian.days <- function(mo, year){
  # detemine how many days are in the month
  if (mo < 12){
    d = as.numeric(difftime(as.Date(paste0(year,"-",as.character(mo+1),"-01")),
                            as.Date(paste0(year,"-",as.character(mo),"-01"))))
  } else {
    d = 31
  }
  # julian days for current month
  day_S = ISOdate(year,mo,1)
  day_E = ISOdate(year,mo,d)
  jd = strptime(seq(day_S,day_E,by='day'), "%Y-%m-%d")$yday+1
  return(jd)
}

# calculate tf models
tf.models <- function(dem0, glacier, months, dir.shift = NULL){
  strt = Sys.time()
  ## Check variables for correct projection
  if(!isLonLat(dem0) | !isLonLat(glacier)){
    stop("projection is not lat/lon")}
  ## Check for similar DEM structure
  try(dem0@extent@xmin, stop("\n !!ERROR: Unknown DEM structure \n  Replace all 'dem@extent@' variables in code"))
  ## VARIABLES (STAY CONSTANT OVER TIME)
  # trim rasters for specific glacier
  dem <- crop.raster(dem0, glacier)
  if(ncell(dem0) == ncell(dem)){
    stop("input DEM extent should be larger (for cast shadow effect)")}
  dmat = as.matrix(dem)
  # slope and aspect (8=queen case, 4=rook case)
  s <- terrain(dem, opt='slope',unit='radians',neighbors=8)
  a <- terrain(dem, opt='aspect',unit='radians',neighbors=8)
  # central coordinates for GLACIER
  lat <- round((dem@extent@ymax + dem@extent@ymin)/2,5); lat_rad = lat * pi/180
  lon <- round((dem@extent@xmax + dem@extent@xmin)/2,5); lon_rad = lon * pi/180
  if(is.na(lat)){stop("Error: DEM has different lat/lon structure")}
  # Z FACTOR and RESOLUTION BASED ON LATITUDE
  Z = (1/(111320*cos(radians(lat))))
  RESOL_DEM <- res(dem)[1]/Z
  # SOLAR CONSTANT (Wm^-2) 
  SC = 1368 
  
  
  ###################################################################################################
  ### MONTH LOOP #####
  # MONTH MATRIX ZEROS
  P1MO = matrix(0,nrow=dim(dmat)[1]*dim(dmat)[2],ncol=length(months))
  P2MO = matrix(0,nrow=dim(dmat)[1]*dim(dmat)[2],ncol=length(months))
  P3MO = matrix(0,nrow=dim(dmat)[1]*dim(dmat)[2],ncol=length(months))
  P4MO = matrix(0,nrow=dim(dmat)[1]*dim(dmat)[2],ncol=length(months))
  for (mo in months){
    print(paste("***beginning month", mo))
    # junlian days for month mo
    jd = julian.days(mo, year)
    # Solar Declination (unit: degrees)
    SD = 23.45*cos(((360*(jd-172))/365)*pi/180)
    SD = SD * pi/180 # convert to radians
    # HOUR ANGLE CALCULATIONS (calculates )
    ### SIMPLE : theta = 90
    HA = acos(- tan(SD)*tan(lat_rad))
    DAY_HOURS = (HA*180/pi/15)*2
    
    # DAY MATRIX ZEROS
    P1D = matrix(0,nrow=dim(dmat)[1]*dim(dmat)[2],ncol=length(jd))
    P2D = matrix(0,nrow=dim(dmat)[1]*dim(dmat)[2],ncol=length(jd))
    P3D = matrix(0,nrow=dim(dmat)[1]*dim(dmat)[2],ncol=length(jd))
    P4D = matrix(0,nrow=dim(dmat)[1]*dim(dmat)[2],ncol=length(jd))
    #################################################################################################
    # DAY LOOP
    for (i in 1:length(jd)){
      print(paste(">>> Starting month",mo,'-',i,"of",length(jd), "days -- at Lat of", lat,"deg"))
      
      # 15 MIN TIME INTERVALS
      TIME_INTERVAL = (15/4)*pi/180 # 4 per hour / 15deg/4 (convert for radians)
      HA_LIST = rev(seq(-HA[i],HA[i], by= TIME_INTERVAL))
      # SOLAR ZENITH & AZIMUTH (IQBAL 2012)
      ZENITH = acos((sin(lat_rad)*sin(SD[i])) + 
                      (cos(lat_rad)*cos(SD[i])*cos(HA_LIST)))*180/pi
      AZIMUTH = acos(((sin((90-ZENITH)*pi/180)*sin(lat_rad)-sin(SD[i])))/
                       (cos((90-ZENITH)*pi/180)*cos(lat_rad)))*180/pi
      # SOLAR NOON
      SOL_NOON = which.min(AZIMUTH)
      # AZIMUTH FOR ZSLOPE CALC
      Az = c(AZIMUTH[1:SOL_NOON]*-1,AZIMUTH[(SOL_NOON+1):length(AZIMUTH)])
      # Change to 360 degrees
      AZIMUTH = c(180 - AZIMUTH[1:SOL_NOON], 180 + AZIMUTH[(SOL_NOON+1):length(AZIMUTH)])
      # CHANGE ASPECT?
      if (is.null(dir.shift)){
        AZIMUTH = change.aspect(AZIMUTH, dir.shift)
      }
      
      #MOMENT MATRIX ZEROS
      P1M = matrix(0,nrow=dim(dmat)[1]*dim(dmat)[2],ncol=length(ZENITH))
      P2M = matrix(0,nrow=dim(dmat)[1]*dim(dmat)[2],ncol=length(ZENITH))
      P3M = matrix(0,nrow=dim(dmat)[1]*dim(dmat)[2],ncol=length(ZENITH))
      P4M = matrix(0,nrow=dim(dmat)[1]*dim(dmat)[2],ncol=length(ZENITH))
      
      # Calculate irradiance at each moment m throughout day i 
      for (m in 1:length(ZENITH)){
        print(paste("making shade for...",m, "of",length(ZENITH),
                    "moments -- Month",mo," Day",i,"- at Lat of",round(lat,2), 'deg'))
        
        # Shading algorythm (Corripio, 2003)
        sv = normalvector(ZENITH[m],AZIMUTH[m])
        sh <- doshade(dmat,sv,dl=RESOL_DEM)
        
        ############################################################
        ## Hock (2005)
        # Earth-sun Radius 
        Rm = 149.6 # semi-major axis (a)
        e = 0.017 # 0.0167
        thet0 = (365.5/360)*(jd[i])
        R = ((Rm*(1-e*e))/(1+e*cos(thet0*pi/180)))
        # Atmospheric Pressure 
        Po = 101325 # atmospheric pressure at sea-level (pa)
        To = 288 # K (15 deg C)
        Lo =-0.0065 # Standard temperature lapse rate (-0.0065 k/m)
        H = 5550 # height about sea level
        Hb = 0 # height at bottom of atmospheric layer
        Rgas = 8.31432 # universal gas constant
        Go = 9.80665 # gravity
        M = 0.0289644 # molar mass of Earth's air
        P = Po*(1+(Lo/To)*(H-Hb))**((-Go*M)/(Rgas*Lo))
        # Atmospheric Transmissivity
        AT = 0.75
        # Equation
        Io = SC*((Rm/R)**2)*AT**(P/(Po*cos(ZENITH[m]*pi/180)))
        ############################################################
        # SOLAR ZENITH AND AZIMUTH ON A SLOPE
        Zenith_rad  = ZENITH[m] * pi/180
        Slope_rad   = s
        Azimuth_rad = Az[m] * pi/180
        Aspect_rad  = a - (180* pi/180)
        THETA = acos((cos(Zenith_rad) * cos(Slope_rad)) +
                       (sin(Zenith_rad) * sin(Slope_rad) * cos(Azimuth_rad - Aspect_rad)))
        # ANGLE OF INCIDENCE MATRIX
        THETA[values(THETA) > 90*pi/180] <- 90*pi/180 
        THETAmat = as.matrix(THETA)
        # ZENITH ANGLE MATRIX
        ZEN = raster(dem)
        ZEN[] <- ZENITH[m]*pi/180
        ZENmat = as.matrix(ZEN)
        
        ###################################################################
        # 4 MODELS 
        P1 = Io*cos(THETAmat)
        P2 = Io*cos(ZENmat)*sh      
        P3 = Io*cos(THETAmat)*sh
        P4 = Io*cos(ZENmat)
        P1M[,m] = array(P1)
        P2M[,m] = array(P2)
        P3M[,m] = array(P3)
        P4M[,m] = array(P4)
        #####################################################################
      }
      # MEAN INSOL FOR DAY
      P1D[,i] = rowMeans(P1M, na.rm = T)
      P2D[,i] = rowMeans(P2M, na.rm = T)
      P3D[,i] = rowMeans(P3M, na.rm = T)
      P4D[,i] = rowMeans(P4M, na.rm = T)
    }
    ########## END DAY LOOP ###########
    # SUM UP MONTHLY VALUES
    P1MO[,match(mo,months)] = rowSums(P1D*is.finite(P1D),na.rm=TRUE)
    P2MO[,match(mo,months)] = rowSums(P2D*is.finite(P2D),na.rm=TRUE)
    P3MO[,match(mo,months)] = rowSums(P3D*is.finite(P3D),na.rm=TRUE)
    P4MO[,match(mo,months)] = rowSums(P4D*is.finite(P4D),na.rm=TRUE)
    # ADD UP TOTAL INSOL HOURS FOR MONTH
    if (mo == months[1]){
      MONTH_HOURS = sum(DAY_HOURS)
      jdays = jd
    } else{
      MONTH_HOURS = c(MONTH_HOURS,sum(DAY_HOURS))
      jdays = c(jdays,jd)
    }
    
  }
  #### END MONTH LOOP #####
  #########################################################################
  ### ADDED MODELS
  # DIVIDE BY NUMBER OF DAYS IN MONTH
  P1MOD = matrix((rowSums(P1MO*is.finite(P1MO),na.rm=TRUE)/length(jdays)),nrow=nrow(dem),ncol=ncol(dem))
  P2MOD = matrix((rowSums(P2MO*is.finite(P2MO),na.rm=TRUE)/length(jdays)),nrow=nrow(dem),ncol=ncol(dem))
  P3MOD = matrix((rowSums(P3MO*is.finite(P3MO),na.rm=TRUE)/length(jdays)),nrow=nrow(dem),ncol=ncol(dem))
  P4MOD = matrix((rowSums(P4MO*is.finite(P4MO),na.rm=TRUE)/length(jdays)),nrow=nrow(dem),ncol=ncol(dem))
  # CREATE STACK TO SAVE
  S_A <- make.raster((P3MOD - P2MOD), dem)
  CS <- make.raster((P3MOD - P1MOD), dem)
  TOT <- make.raster((P3MOD - P4MOD), dem)
  SR <- make.raster(((P2MOD - P4MOD) - (P3MOD - P1MOD)), dem) # double check
  ###########################################################################
  # SAVE ALL VALUES TO STACK
  tmp_stk = stack(S_A, CS, SR, TOT, dem)
  names(tmp_stk) = c("Slope and Aspect", "Cast Shadow", "Shaded Relief", "Combined", "DEM")
  
  # redefine extent for better images
  # Create Extent
  ex = extent(glacier)
  ex@ymin = ex@ymin-0.001
  ex@ymax = ex@ymax+0.001
  tf_models = mask(crop(tmp_stk, ex), glacier)
  # save
  if (!is.null(save_path)){
    spath = paste0(save_path,"/TopoSol")
    dir.create(spath)
    # writeRaster(tf_models, filename = paste0(spath,"/tf_models.grd"), format="raster")           # R file
    writeRaster (tf_models, filename = paste0(spath,"/tf_models.tif"), options = c('TFW = YES'))   # ESRI GeoTiff
  }
  ##########################################################################
  #### END LAT LOOP #####
  end = Sys.time() - strt 
  print(end)
  return(tf_models)
}


