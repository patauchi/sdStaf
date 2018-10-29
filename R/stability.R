#' Stability of ecological niche models

stability <- function(current=NULL, project=NULL, thr.value,
                      continue = FALSE)
#'
#' Returns stability niche based on species distribution models and their projections.
#'
#'
#' @param current Raster* objet of present distribution.
#' Raster has continue values.
#' @param project RasterStack* object of project distributions.
#' Must have three models with continue values.
#' @param thr.value Cut value (0 — 1) of threshold in order to
#'  species distribution.
#'
#' @param continue defines if the species distribution are either binary or continue maps.
#' \code{ Default = FALSE}.
#'
#' @return An object of class 'StabEcodist'
#'
#' \bold{Based on binary maps}
#'
#' Return table with these features: \code{Models} and \code{nPixel}
#' (frequency of pixel with that feature).
#'
#' Stability maps based on  binary species distirbution, give:
#'
#'  \code{Values of 0 } Shows absence
#'
#'
#'  \code{Values of 100 } Mentions the lost area
#'
#'
#'  \code{Values [1:100]} Shows colonizable area. Different models are defined as numbers (e.g.
#' Value of \code{1} indicates one models predict gain; Value of \code{2} indicates
#'  two models predict agains)
#'
#'
#'   \code{Values  > 100 } Shows stability or permanence. Differente models are defined as numbers (e.g.
#' value of \code{101} mentions one model predict stability, value of \code{102} mentions
#' two models predict stability)
#'
#'
#'
#' \bold{Based on continue maps}
#'
#' Species distribution show different values of stability along
#' of their distributions.
#'
#'\code{Values of -2 } Shows absence
#'
#'
#'  \code{Values [-1 : 0] } Shows colonized grade or gain.
#'
#'
#'  \code{Values [0 : 1] } Shows stability or permanence
#'
#'
#'    \code{Values of 2 } Shows lost area
#'
#'
#' @references
#' Peterson et al., (2017) Influences of climate change on the potential
#' distribution of Lutzomyia longipalpis sensu lato (Psychodidae: Phlebotominae).
#'  International Journal for Parasitology. 47(10–11):667–74.
#'
#' @export


{

  if(is.null(current) == TRUE){
    stop('You need to define current distribution')
  }
  if(is.null(list(project)) == TRUE){
    stop('You need to define projection distribution')
  }

  if(continue == FALSE){

    if ((thr.value >= 0 & thr.value <= 1) == FALSE) {
      message('Threshold is value between [0 - 1]')
    }
    # Build stability analysis
    dfa <- stabl(current, project, thr.value)

    # data.frame of stability
    stab.an <- raster::rasterToPoints(dfa)
    stab.an <- as.data.frame(stab.an)

    # show statistics
    stab.an <- as.data.frame(table(stab.an$layer))
    names(stab.an) <- c("Models", "nPixel")
    stin <- stab.an

    output <- new("StabEcodist", df = stin, map = dfa)
    return(output)

  }else{

    dfa <- stabl.con(current, project, thr.value)

    output <- new("StabEcodist", map = dfa)
    return(output)
  }

}
