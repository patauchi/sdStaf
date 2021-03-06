
#' Reduce environmental data
#' @description
#' This function allows reduce environmental data clipping by mask or buffer area
#'
#' @details
#' Reduce the correlation among predicted variables either buffer zone, or
#' clipping mask.
#'
#' Provide reduce objet with cut predicted variables and data.frame for
#' the values of each point of occurrence of them.
#'
#' @param env RasterStack* objet.
#' @param transfer List of RasterStack object
#' @param var.rm A data.frame of occurrence records.
#'  It must include two column based on latitude and longitude.
#' @param mask Cropped mask, must be shapefile (.shp), readOGR.
#'
#' @return
#'
#'  \code{@crop} RasterStack* Objet
#'
#'  \code{@m.env} data.frame of environmental values to occurrence localities.
#'
#' @seealso \code{\link{cor.show}}
#'
#' @examples
#'
#' # Phytotoma ocurrence data
#' data(phytotoma)
#'
#'
#' # Complement
#' library(dismo)
#' predictor <- stack(list.files(path=paste(system.file(package="dismo"),'/ex', sep=''),
#'  pattern='grd', full.names=TRUE ))
#'
#'  maskM <- stim.M(phytotoma[,2:3], 131)
#'
#' vrn <- extract(predictor, phytotoma[,2:3])
#'
#' # Correlogram s
#' rd <- c('bio1','bio12','biome')
#' cor.show(vrn, remove=TRUE, var.rm=rd)
#' 
#'
#'
#' @export
#'
reduce.env <- function(env, transfer=NULL, var.rm, mask)
{
  ptm <- proc.time()
  if(is.null(env)){
    stop('You need to define environmental data')
  }
  if (is.null(var.rm)){
    stop('You need to define variables to drop')
  }
  if(is.null(mask)){
    stop('You need to build M hypothesis before to do this')
  }
  
  mask.M <- mask
  
  if(is.null(transfer)){
    
    # DropLayer
    pr_var <- dropLayer(env, var.rm)
    
    #Mask of datalayer
    biovars.mask <- raster::crop(pr_var, mask.M)
    biovars.mask <- raster::mask(biovars.mask, mask.M)
    
    layer.transfer <- list()
    datavalue <- matrix()
  }else{
    # Subset
    pr_var <- raster::subset(env, var.rm)
    fut_var <- list()
    for (i in 1:length(transfer)) {
      fut_var[[i]] <- raster::subset(transfer[[i]], var.rm)
    }
    
    # Mask of datalayer
    biovars.mask <- raster::crop(pr_var, mask.M)
    biovars.mask <- raster::mask(biovars.mask, mask.M)
    
    # Mask of future datalayer
    layer.transfer <- list()
    # Core function
    for (i in 1:length(fut_var)) {
      layer.transfer[[i]] <- raster::crop(fut_var[[i]], mask.M)
      layer.transfer[[i]] <- raster::mask(layer.transfer[[i]], mask.M)
    }
    datavalue <- matrix()
  }
  
  r <- EnvimRed(cropa = biovars.mask,
                 project = layer.transfer,
                 m.env = datavalue)
  
  timed <- proc.time() - ptm
  t.min <- floor(timed[3] / 60)
  t.sec <- timed[3] - (t.min * 60)
  
  message(paste("Environmental data completed in", t.min, "minutes", round(t.sec, 2), "seconds."))

  return(r)

}
