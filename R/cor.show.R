#' Correlogram plot
#'
#' Returns plot with correlation values among predicted variables.
#'
cor.show <- function(values, rm=FALSE, var.rm)

#' Correlation matrix based on pearson.
#'
#' @param values DataFrame-class.
#' @param rm logical. If \code{TRUE}, allows remove some
#' variables from imput data set. (\code{default = FALSE})
#' @param var.rm variables names of imput data set. Using
#'  \code{colnames(RasterStick)}, Where \code{RasterStick} is RasterStack* object.
#'
#' @return Correlogram plot
#'
#' @importFrom raster select
#' @importFrom dplyr one_of
#' @importFrom graphics panel.smooth
#'
#' @seealso \code{\link{reduce.env}}
#'
#' @export
#'
{ # r es la clase producida por la funcion red.env.
   r <- values
  if(rm==FALSE){
    datavalue <- r

  } else{
    #datavalue <- datavalue[ ,!colnames(datavalue) == c('bio_12','bio_19','bio_18','bio_8') ]
    datavalue <- (data.frame(r)) %>% dplyr::select(-one_of(var.rm))
    datavalue <- data.matrix(datavalue, rownames.force = NA)

  }

  # Extrae la matriz con la tabla de datos
  nm <- colnames(datavalue)
  nn <- length(nm)

  # produce el plot de correlaciones y valores de r
  pairs(datavalue[,1 : nn],lower.panel=panel.smooth,
        upper.panel=panel.r2,diag.panel=panel.hist)

}
