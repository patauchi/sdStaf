if (getRversion() >= "2.15.1") { utils::globalVariables(c("sd", "IQR"))}
#' Build buffer zone to M
#' 
#' Returns buffer zone based on ocurrence data


stim.M <- function (occs, radio=NULL, bgeo=NULL, method='user', env=NULL, Vrc = 1, ncal = 1, ...)


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
#' @param env if True. Environmental daataset used to build M. Only \code{method = 'Tol.pca'}
#' @param Vrc Integer. sd(IQR) * value, used to increase range tolerance of dataset \code{env}
#' @param ncal Integer. Dataset using to define IQR. Only \code{method = 'Tol.pca'}
#' @param method default = 'user'. Another option is calculate the mean of all points 'mean'.
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
#' 

{
  
  if(is.null(bgeo)){
    
    METHODS <- c("user", "Mx.dist","mean","Tol.pca")
    
    method <- match.arg(method, METHODS)
    
    if(method == "user")
    {
      if(is.null(radio)){
        stop('Define calibration radio')
      }else{
        
        rat <- 1000 * radio
        sp_po <- SpatialPoints(occs)
        projection(sp_po) <- CRS('+proj=longlat +datum=WGS84')
        hM.pol <- raster::buffer(sp_po, width = rat)
      }
    }
    
    if(method == "Mx.dist"){
      dm <- geosphere::distm(occs)
      radio <- max(dm)/1000
      rat <- 1000 * radio
      sp_po <- SpatialPoints(occs)
      projection(sp_po) <- CRS('+proj=longlat +datum=WGS84')
      hM.pol <- raster::buffer(sp_po, width = rat)
    }
    if(method == "mean"){
      dm <- geosphere::distm(occs)
      radio <- mean(dm)/1000
      rat <- 1000 * radio
      sp_po <- SpatialPoints(occs)
      projection(sp_po) <- CRS('+proj=longlat +datum=WGS84')
      hM.pol <- raster::buffer(sp_po, width = rat)
    }
    
    if(method == "Tol.pca"){
      if(is.null(env)) {stop("Define env argument to build Tolerance Range")}
      if(nlayers(env) <= ncal)
        stop("Dataset define to build M is so longer")
      
      if(is.null(radio)){
        dm <- geosphere::distm(occs)
        radio <- mean(dm)/1000
        rat <- 1000 * radio
        sp_po <- SpatialPoints(occs)
        projection(sp_po) <- CRS('+proj=longlat +datum=WGS84')
        spbuf <- raster::buffer(sp_po, width = rat)
        
        c1 <- raster::crop(env, spbuf)
        c1 <- raster::mask(c1, spbuf)
        
      }else{
        rat <- 1000 * radio
        sp_po <- SpatialPoints(occs)
        projection(sp_po) <- CRS('+proj=longlat +datum=WGS84')
        spbuf <- raster::buffer(sp_po, width = rat)
        
        c1 <- raster::crop(env, spbuf)
        c1 <- raster::mask(c1, spbuf)
        
      }
      
      gt_p <- as.data.frame(extract(c1, occs))
      
      contn <- list()
      
      for (i in 1:nlayers(c1)) {
        
        LsupC = quantile(gt_p[[i]], probs = 0.75) + 1.5 * (IQR(gt_p[[i]])) + sd(gt_p[[i]]) * Vrc
        LinfC =  quantile(gt_p[[i]], probs = 0.25) - 1.5 * (IQR(gt_p[[i]])) - sd(gt_p[[i]]) * Vrc
        
        s1 <- ifelse(getValues(c1[[i]]) > LinfC & getValues(c1[[i]]) < LsupC, 1, 0 )
        m1 <- setValues(c1[[i]], s1)
        
        contn[[i]] <- m1 
      }; contn <- stack(contn)
      
      contn <- calc(contn, sum)
      #plot(contn)
      
      s1 <- ifelse(getValues(contn) >= ncal, 1, 0)
      # Re
      hM.ras <- setValues(contn, s1)
      #plot(hM)
      
      hM.pol <- rasterToPolygons(hM.ras, fun=function(x){x==1}, dissolve = T) 
      #plot(hM.pol)
      #return(hM.pol)
    }
    return(hM.pol)
  }else{
    
    b_ex <- raster::extract(bgeo, occs)
    b_ex <- unique(b_ex$ECO_NAME)
    b_ex <- as.vector(b_ex)
    
    shapeOut <- subset(bgeo, bgeo@data[,4] %in% b_ex)
    projection(shapeOut) <- CRS('+proj=longlat +datum=WGS84')
    #rat <- 1000 * radio
    #spbuf <- buffer(sp_po, width = rat)
    
    
    METHODS <- c("user", "Mx.dist","mean","Tol.pca")
    
    method <- match.arg(method, METHODS)
    
    if(method == "user")
    {
      if(is.null(radio))
        stop('Define calibration radio')
      radio
      
      rat <- 1000 * radio
      sp_po <- SpatialPoints(occs)
      projection(sp_po) <- CRS('+proj=longlat +datum=WGS84')
      spbuf <- raster::buffer(sp_po, width = rat)
      
      hM.pol <- raster::intersect(shapeOut, spbuf)
      
    }
    
    if(method == "Mx.dist"){
      dm <- geosphere::distm(occs)
      radio <- max(dm)/1000
      
      rat <- 1000 * radio
      sp_po <- SpatialPoints(occs)
      projection(sp_po) <- CRS('+proj=longlat +datum=WGS84')
      spbuf <- raster::buffer(sp_po, width = rat)
      
      hM.pol <- raster::intersect(shapeOut, spbuf)
    }
    if(method == "mean"){
      dm <- geosphere::distm(occs)
      radio <- mean(dm)/1000
      
      rat <- 1000 * radio
      sp_po <- SpatialPoints(occs)
      projection(sp_po) <- CRS('+proj=longlat +datum=WGS84')
      spbuf <- raster::buffer(sp_po, width = rat)
      
      hM.pol <- raster::intersect(shapeOut, spbuf)
    }
    
    
    if(method == "Tol.pca"){
      if(is.null(env))
        stop("Define env argument to build Tolerance Range")
      if(nlayers(env) <= ncal)
        stop("Dataset define to build M is so longer")
      
      if(is.null(radio)){
        c1 <- raster::crop(env, shapeOut)
        c1 <- raster::mask(c1, shapeOut)
        
      }else{
        rat <- 1000 * radio
        sp_po <- SpatialPoints(occs)
        projection(sp_po) <- CRS('+proj=longlat +datum=WGS84')
        projection(shapeOut) <- CRS('+proj=longlat +datum=WGS84')
        spbuf <- raster::buffer(sp_po, width = rat)
        
        temp.M <- raster::intersect(shapeOut, spbuf)
        
        c1 <- raster::crop(env, temp.M)
        c1 <- raster::mask(c1, temp.M)
        
      }
      
      gt_p <- as.data.frame(extract(c1, occs)) %>% na.omit()
      
      contn <- list()
      
      for (i in 1:nlayers(c1)) {
        
        LsupC = quantile(gt_p[[i]], probs = 0.75) + 1.5 * (IQR(gt_p[[i]])) + sd(gt_p[[i]]) * Vrc
        LinfC =  quantile(gt_p[[i]], probs = 0.25) - 1.5 * (IQR(gt_p[[i]])) - sd(gt_p[[i]]) * Vrc
        
        s1 <- ifelse(getValues(c1[[i]]) > LinfC & getValues(c1[[i]]) < LsupC, 1, 0 )
        m1 <- setValues(c1[[i]], s1)
        
        contn[[i]] <- m1 
      }; contn <- stack(contn)
      
      
      contn <- calc(contn, sum)
      #plot(contn)
      
      s1 <- ifelse(getValues(contn) >= ncal, 1, 0)
      # Re
      hM.ras <- setValues(contn, s1)
      #plot(hM)
      
      hM.pol <- rasterToPolygons(hM.ras, fun=function(x){x==1}, dissolve = T) 
      #plot(hM.pol)
      
      
      #message(paste0('We have defined the M radio based on', method , 'method'))
      
      
    }
  } 
  return(hM.pol)
}


