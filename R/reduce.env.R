
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
#'
reduce.env <- function(env, transfer=NULL, occ_data, mask)
{
  ptm <- proc.time()
  if(is.null(env)){
    stop('You need to define environmental data')
  }
  if (is.null(occ_data)){
    stop('You need to define ocurrence data (Longitude/Latitude)')
  }

  if(is.null(transfer)){
    # Corta las variables ambientales originales (todo el mundo) hacia el area de interes.
    biovars.mask <- crop (env, mask)
    biovars.mask <- mask (biovars.mask, mask)
    layer.transfer <- list('Don not have environment data')
    datavalue <- extract(env, occ_data)
    datavalue <- na.omit(datavalue)

  } else{

    biovars.mask <- crop (env, mask)
    biovars.mask <- mask (biovars.mask, mask)
    layer.transfer <- list()
    for (i in 1:length(transfer)) {
      layer.transfer[[i]] <- crop(transfer[[i]], mask)
      layer.transfer[[i]] <- mask(layer.transfer[[i]], mask)
    }
    datavalue <- extract(biovars.mask, occ_data)
    datavalue <- na.omit(datavalue)

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
