if (getRversion() >= "2.15.1") { utils::globalVariables(c("detectCores","makeCluster",
                                                          "registerDoParallel","%dopar%","foreach",
                                                          "stopCluster","distm"))}

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
#' @param transfer List of rasterstack object
#' @param occ_data A data.frame of occurrence records.
#'  It must include two column based on latitude and longitude.
#' @param mask Croped mask, must be shapefile (.shp), readOGR.
#' @param parallel Logical. Build parallel process with each future project. Default is FALSE.
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
#' reduce_cut <- reduce.env(env = predictor, occ_data = phytotoma[,2:3], mask=maskM)
#'
#' # Plot reduce_cut
#' plot(reduce_cut@cropa$bio1)
#'
#' # Add points
#' points(phytotoma[,2:3], pch=16,col='blue')
#'
#' # Correlogram
#' cor.show(reduce_cut)
#' rd <- c('bio1','bio12','bio16','biome','bio8')
#'
#' # Removing rd-variables on correlogram
#' cor.show(reduce_cut, rm=TRUE, var.rm = rd)
#'
#' # Remove rd-variables
#' var_reduce <- dropLayer(reduce_cut@cropa, rd)
#'
#' # summary
#' var_reduce
#'
#'
#' @export
#'
reduce.env <- function(env, transfer=NULL, occ_data, mask, parallel = FALSE)
{
  ptm <- proc.time()
  if(is.null(env)){
    stop('You need to define environmental data')
  }
  if (is.null(occ_data)){
    stop('You need to define ocurrence data (Longitude/Latitude)')
  }
  if(is.null(mask)){
    stop('You need to build M hypothesis before to do this')
  }

  if(is.null(transfer)){
    if(parallel == TRUE){
      message('You have a one dataset, parallel = FALSE ' )
    }
    mask.M <- mask
    # Corta las variables ambientales originales (todo el mundo) hacia el area de interes.
    biovars.mask <- crop (env, mask.M)
    biovars.mask <- mask (biovars.mask, mask.M)
    layer.transfer <- list('Don not have environment data')
    datavalue <- extract(env, occ_data)
    datavalue <- na.omit(datavalue)

  } else{
    if(parallel == FALSE){
      biovars.mask <- crop (env, mask.M)
      biovars.mask <- mask (biovars.mask, mask.M)
      
      layer.transfer <- list()
      
      # Core function
      for (i in 1:length(transfer)) {
        layer.transfer[[i]] <- crop(transfer[[i]], mask.M)
        layer.transfer[[i]] <- mask(layer.transfer[[i]], mask.M)
      }
      # Timer produce
      # tock    <- proc.time()[3]
      #  timeG[1] <- tock - tick
      # --- End Timer
      
      datavalue <- extract(biovars.mask, occ_data)
      datavalue <- na.omit(datavalue)  
    } else {
      biovars.mask <- crop (env, mask.M)
      biovars.mask <- mask (biovars.mask, mask.M)
      
      #library(foreach)
      #library(doParallel)
      
      cores = detectCores()
      cl <- makeCluster(cores[1]-1) #not to overload your computer
      
      registerDoParallel(cl)
      
      layer.transfer <- list()

      # Core function
      layer.transfer <- foreach(i=1:length(transfer), .packages = 'raster') %dopar% {
        
        layer.transfer[[i]] <- crop(transfer[[i]], mask.M)
        layer.transfer[[i]] <- mask(layer.transfer[[i]], mask.M)
      };stopCluster(cl)
      
      
      datavalue <- extract(biovars.mask, occ_data)
      datavalue <- na.omit(datavalue)
    }
    

  }

  
  r <- EnvimRed(cropa = biovars.mask,
                 project = layer.transfer,
                 m.env = datavalue)

  timed <- proc.time() - ptm
  t.min <- floor(timed[3] / 60)
  t.sec <- timed[3] - (t.min * 60)
  
  paste("Environmental data completed in", t.min, "minutes", round(t.sec, 2), "seconds.")

  return(r)

}
