if (getRversion() >= "2.15.1") { utils::globalVariables(c("coordinates<-", "proj4string<-","over","mutate",
                                                          "Specie","Longitude","Latitude","fac"))}

#' Thin Spatial Occurrence
#'
#' Returns a data frame of occurrence dataset.
#'
ThinOcc <- function(df, NullRaster, SpeciesName=NULL, ThinRange=1, Polygon=NULL,
                    methods=c('spatial','polygons'), save=TRUE,...)
  
#' Thin ocurrence data using three approaches Polygons, Buffer zone, and Biogeographical regions.
#'
#' @param df dataframe-class. Occurrence data e.g. long,lat.
#' @param NullRaster raster-class. 
#' @param SpeciesName character-class.
#' @param ThinRange numeric-class. Radio of filter species occurrence.
#' @param Polygon spatial-class.
#' @param methods Choose .
#' @param save logical. 
#' @param ... no implemented. Use \code{trace=TRUE} to show vif running. \code{default=FALSE}.
#'
#' @seealso \code{\link{reduce.env}}
#' 
#' @importFrom stats complete.cases
#' @importFrom utils write.csv
#' 
#' @export
{
  Occ <- df # <- dfXY
  env <- NullRaster#  <- predictor$bio1 
  sp.name <- SpeciesName  # <- 'sp1'
  
  if(methods=='polygons'){
    if(is.null(Polygon)){stop(message('Polygon parameter is empty. Make sure that you are using a correct parameter.'))}
    
    colnames(Occ) <- c('X','Y')
    coordinates(Occ) <- ~X+Y
    sps <- Polygon#  <- sps
    
    proj4string(sps) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    proj4string(Occ) <-  CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    
    dfOccB <- Occ[complete.cases(over(Occ, sps)), ]
  }
  
  if(methods=='spatial'){
    dfOccB <- Occ

  }

  
  # Subset
  sp.d <- sp::SpatialPoints(dfOccB)
  raster::crs(sp.d) <- c('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
  
  sp.d <- sp::coordinates(sp.d)
  
  # Ajuste de ocurrencias
  data_occ <- spThin::thin.algorithm(sp.d, thin.par = ThinRange, reps = 1)
  data_occ <- data_occ[[1]]
  
  ## crear extent
  lco <- (sp::SpatialPoints(data_occ))
  
  coords = matrix(c(raster::extent(lco)[1] - 0.5, raster::extent(lco)[3] - 0.5,
                    raster::extent(lco)[2] + 0.5, raster::extent(lco)[3] - 0.5,
                    raster::extent(lco)[2] + 0.5, raster::extent(lco)[4] + 0.5,
                    raster::extent(lco)[1] - 0.5, raster::extent(lco)[4] + 0.5), 
                  ncol = 2, byrow = TRUE)
  
  
  pM = sp::Polygon(coords)
  
  pM = sp::SpatialPolygons(list(sp::Polygons(list(pM), ID = "a")), 
                           proj4string=sp::CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
  #pM = SpatialPolygonsDataFrame(list(Polygons(list(pM), ID = "a")))
  #plot(pM)
  
  env <- raster::crop(env,pM)
  env <- raster::mask(env,pM)
  raster::crs(env) <- c('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
  
  
  ### Create background points
 suppressWarnings(bg.pts <- dismo::randomPoints(env, 500))
  
  ### Checkerboard1 partitioning method
  chk1.pts <- ENMeval::get.randomkfold(data_occ, bg.pts, 3)
  
  #points(data_occ, pch=23, bg=chk1.pts$occ.grp)
  
  #length(chk1.pts$occ.grp[chk1.pts$occ.grp==1])
  #length(chk1.pts$occ.grp[chk1.pts$occ.grp==2])
  #length(chk1.pts$occ.grp[chk1.pts$occ.grp==3])
  
  data_occ2 <- data.frame(data_occ, fac = chk1.pts$occ.grp)#
  
  train <- data_occ2 %>%
    dplyr::filter(fac != 1) %>%
    dplyr::mutate(Specie = sp.name) %>%
    dplyr::select(Specie, Longitude, Latitude) %>%
    data.frame()
  
  test <- data_occ2 %>%
    dplyr::filter(fac == 1) %>%
    dplyr::mutate(Specie = sp.name) %>%
    dplyr::select(Specie, Longitude, Latitude) 
  
  joint <- data_occ2 %>%
    dplyr::mutate(Specie = sp.name) %>%
    dplyr::select(Specie, Longitude, Latitude) 
  
  
  if(save){
    write.csv(test, 'sp_test.csv', row.names = F)
    write.csv(train, 'sp_train.csv', row.names = F)
    write.csv(joint, 'sp_joint.csv', row.names = F)
    
  } else if(save==FALSE){
    
    Output <- list()
    Output$Joint <- joint
    Output$Train <- train
    Output$Test <- test
    
    return(Output)
    
  }
  
  message('We have ended create occurrence dataset for ', sp.name) 
  
}
