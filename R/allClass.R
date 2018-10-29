if (getRversion() >= "2.15.1") { utils::globalVariables(c("value"))}
#'  myRefClass

#' @import raster dismo methods ggplot2
#' @importFrom raster select union stack density hist lines pairs text
#' @importFrom dplyr %>%
#' @importFrom ggplot2 geom_tile aes scale_fill_manual ggplot
#' @importFrom rasterVis gplot
#'
#' @exportClass EnvimRed StabEcodist
#'
#' EnvimRed
#' @name EnvimRed-class
#' @rdname EnvimRed-class
#' @slot cropa A RasterBrinck
#' @slot m.env A Matrix
#' @slot project A List
EnvimRed <- setClass("EnvimRed",
                      slots = c(cropa="RasterBrick",
                                project = "list",
                                m.env="matrix"))

#'  StabEcodist
#' @name StabEcodist-class
#' @rdname StabEcodist-class
#' @slot df A RasterBrinck
#' @slot map A Matrix
#'

setClass("StabEcodist",
                      slots = c(df="data.frame",
                                map="RasterLayer"))


setMethod("print","StabEcodist",
           function(x, ...){
            cat("*** Class Trajectories, method Print *** \n")
            cat("* Times =");
            print (x@df)
            cat("******* End Print (trajectories) ******* \n")
          })




setMethod("plot","StabEcodist",
          definition = function(x, y="missing", ...){

            mape <- x@map
            rasterVis::gplot(mape) +
              ggplot2::geom_tile(aes(fill = factor(value, labels = c("Absence","Gain 1 model", "Gain 2 model",
                                     "Gain 3 model", "Lost", "Stab 1 model",
                                     "Stab 2 model", "Stab 3 model"))),
                alpha=0.8)+
              scale_fill_manual(values = c("gray", "yellow2","orange2","red",
                                           "limegreen","skyblue","dodgerblue2",
                                           "blue"),
                                name= "Stability map")
          })

requireNamespace("sp", quietly=TRUE)
requireNamespace("raster", quietly=TRUE)
requireNamespace("rJava", quietly=TRUE)
