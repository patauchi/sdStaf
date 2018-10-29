#' Build buffer zone to M
#'
#' Returns buffer zone based on ocurrence data

stim.M <- function (occs, radio, bgeo=NULL, ...)

#'
#' To define calibration area is crucial step (Barve et al., 2011),
#' even more with incomplete sample data sometime is
#' complicated, because to get complete sample within geography space
#' is dificult, in these cases is apropiate define M with buffer zone
#'  (Peterson et al., 2017); and in other cases it helps to cut the
#' ends of the calibration area based on the maximum dispersion capacity
#' (Atauchi et al., 2018).
#'
#' @param occs data.frame of ocurrence data (longitude/latitude).
#' @param radio radio of buffer.
#' @param bgeo Biogeographical layer. Categorical values.
#' @param ... Optional features of buffer
#'
#' @return SpatialPolygons* object
#'
#'
#' @references
#' Atauchi et al. (2018). Species distribution models for
#' Peruvian Plantcutter improve with consideration of biotic
#' interactions. \emph{J. avian biology 2018: e01617}. <doi:http://10.1111/jav.01617.>
#'
#' Barve et al. (2011) The crucial role of the accessible area in
#' ecological niche modeling and species distribution modeling.
#'  \emph{Ecol. Mod}. 222:1810–1819.
#'
#' Peterson et al.(2017) Influences of climate change on the potential
#' distribution of \emph{Lutzomyia longipalpis} sensu lato (Psychodidae:
#' Phlebotominae). \emph{International journal for parasitology}.
#' 45(10-11): 667–674.
#'
#' @importFrom sp CRS SpatialPoints
#'
#' @examples
#'
#' # Phytotoma ocurrence data
#' data(phytotoma)
#'
#' # Build buffer zone
#' buf_M <- stim.M(occs=phytotoma[,2:3], 100)
#'
#' # Add points
#' points(phytotoma[,2:3])
#'
#' @export
#'
{
  if(is.null(bgeo)){
    rat <- 1000 * radio
    sp_po <- SpatialPoints(occs)
    projection(sp_po) <- CRS('+proj=longlat +datum=WGS84')
    spbuf <- buffer(sp_po, width = rat)
    plot(spbuf)
    return(spbuf)


  }else{

    b_ex <- raster::extract(bgeo, occs)
    b_ex <- unique(b_ex$Province_1)
    b_ex <- as.vector(b_ex)

    shapeOut <- subset(bgeo, bgeo@data[,4] %in% b_ex)

    #plot(shapeOut)
    #points(occ[,2:3], pch=16,cex=0.5)
    rat <- 1000 * radio
    sp_po <- SpatialPoints(occs)
    projection(sp_po) <- CRS('+proj=longlat +datum=WGS84')
    spbuf <- buffer(sp_po, width = rat)

    #plot(bu_z)
    #points(occ[,2:3], pch=16,cex=0.5)

    spbuf <- crop(shapeOut, spbuf)
    plot(spbuf)
    #plot(M.zo)
    #paste0('Mask based on Biogeography. Morrone (2014)')
    #paste0(b_ex)
    return(spbuf)
  }








}
