
#' @importFrom graphics par strwidth
#' @importFrom stats cor na.omit reorder filter

##################  AUTO-CORRELATION   #############
# Create function to make histogram with density superimposed
panel.hist <- function(x, ...) {
  # Set user coordinates of plotting region
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  par(new=TRUE)
  # Do not start new plot
  hist(x, prob=TRUE, axes=FALSE, xlab="", ylab="",
       main="", col="lightgray")
  lines(density(x, na.rm=TRUE))
  # Add density curve
}
# Create function to compute and print R^2
panel.r2 <- function(x, y, digits=2, cex.cor, ...) {
  # Set user coordinates of plotting region
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use="complete.obs")**2 # Compute R^2
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  if(missing(cex.cor)) cex.cor <- 1/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}


stabl <- function(p, q, thr.value)

{
  # build threshould
  t <- as.numeric(thr.value)
  t <- c(0, t, 0, t, 1, 1)
  pm <- c(0, 0, 0,  0, 1, 100)

  # build matrix
  mr <- matrix(t, ncol = 3, byrow = TRUE)
  rcls <- matrix(pm, ncol=3, byrow=TRUE)
  p <- reclassify(p, mr)
  pmf <- reclassify(p, rcls)

  qs <- reclassify(q, mr, right=FALSE)
  qs <- calc(qs, fun=sum)

  # Build stability
  stab <- pmf + qs
  #mapp <- rasterToPoints(stab)

  #df <- data.frame(mapp)
  #colnames(df) <- c("Longitude", "Latitude", "MODEL")
  return(stab)
}


stabl.con <- function(p, q, thr.value){
  # build threshould
  t <- as.numeric(thr.value)
  ts <- c(0, t, -1)
  mr <- matrix(ts, ncol = 3, byrow = TRUE)
  c_pre <- reclassify(p, mr,right=FALSE)

  t1 <- c(0, t, 0)
  mr1 <- matrix(t1, ncol = 3, byrow = TRUE)
  qs <- reclassify(q, mr1,right=FALSE)

  fun <- function(x){
      # paramaters
        nparm <- length(x[x > 0])
      # core function
        mean(x) * nparm / nlayers(q) }

  q.val <- calc(qs, fun = fun)
    v <- getValues(q.val)
    v <- replace(v, v == 0, 10)

    qs.complex <- setValues(q.val, v)

  plot(qs.complex)


  # Build stability
  stab <- c_pre * qs.complex

  t3 <- c(-10, -1, -2, 1 , 10, 2)
  mr3 <- matrix(t3, ncol = 3, byrow = TRUE)
  stab.sal <- reclassify(stab, mr3, right=FALSE)


  #plot(stab.sal)
  #con.rast <- stack(q.val, c_pre)
  #mX.lost <- function(x, y, q = 2, p = -1) {ifelse(x < 0.1, p * q, x * y)}
  #stab <- overlay(x=q.val,y=c_pre, fun = mX.lost, unstack=TRUE,forcefun=T)

#writeRaster(stab.sal, 'Prueb.asc',overwrite=T)

  #plot(stab)
  #df <- data.frame(mapp)
  #colnames(df) <- c("Longitude", "Latitude", "MODEL")
  return(stab.sal)
}



