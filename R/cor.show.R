#' Correlogram plot
#'
#' Returns plot with correlation values among predicted variables.
#'
cor.show <- function(values, remove=FALSE, var.rm)

#'  Matrix based on Pearson correlation.
#'
#' @param values DataFrame-class.
#' @param remove logical. If \code{TRUE}, allows remove some
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
#' @examples
#' 
#' DFset <- data.frame(
#' var1 = c(1,2,3,2,3,2,3,2,3,2,3,4,3,4,4),
#' var2 = c(2,3,6,4,5,4,5,6,4,2,3,4,3,6,1),
#' var3 = c(5,4,3,4,5,4,5,6,5,6,2,3,4,5,3),
#' var4 = c(5,4,3,1,2,3,5,3,3,3,2,4,3,3,4),
#' var5 = c(5,2,5,1,1,2,4,1,2,3,3,4,1,2,5),
#' var6 = c(10,12,23,34,23,34,23,12,9,23,12,34,12,23,9))
#' 
#' cor.show(DFset)
#' 
#'
#' @export
#'
{ # r es la clase producida por la funcion red.env.
   r <- values
  if(remove==FALSE){
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
        upper.panel=panel.r2,diag.panel=panel.hist,pch=16)

}
