
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

#  plot(qs.complex)


  # Build stability
  stab <- c_pre * qs.complex

  t3 <- c(-10, -1, -2, 1 , 10, 2)
  mr3 <- matrix(t3, ncol = 3, byrow = TRUE)
  stab.sal <- reclassify(stab, mr3, right=FALSE)

  return(stab.sal)
}


require(clusterGeneration)
VIF_Selection<-function(df.Stack,thresh,trace=F,...){
  
  in_frame <- df.Stack
  #library(fmsb)
  if(any(!'data.frame' %in% class(in_frame))) in_frame<-data.frame(in_frame)
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  var_names <- names(in_frame)
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- stats::formula(paste(val, '~', form))
    vif_init<-rbind(vif_init, c(val, fmsb::VIF(stats::lm(form_in, data = in_frame, ...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]), na.rm = TRUE)
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(var_names)
  }
  else{
    in_dat<-in_frame
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      vif_vals<-NULL
      var_names <- names(in_dat)
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in <- stats::formula(paste(val, '~', form))
        vif_add<-fmsb::VIF(stats::lm(form_in, data = in_dat, ...))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2]), na.rm = TRUE))[1]
      vif_max<-as.numeric(vif_vals[max_row,2])
      if(vif_max<thresh) break
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        utils::flush.console()
      }
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
    }
    return(names(in_dat))
  }
}

### CORRELATION
MTHs_Heuristico <- function (df.Stack, multicollinearity.cutoff, 
                             select.variables = TRUE, sample.points = FALSE, 
                             plot = FALSE, method = "pearson") {
  
  if (!is.numeric(multicollinearity.cutoff)) {
    stop("You must provide a numeric cutoff between 0 and 1 in multicollinearity.cutoff")
  }
  else if (multicollinearity.cutoff > 1 | multicollinearity.cutoff < 0) {
    stop("You must provide a numeric cutoff between 0 and 1 in multicollinearity.cutoff")
  }
  #cor.matrix <- matrix(data = 0, nrow = nrow(df.Stack), 
   #                    ncol = ncol(df.Stack), dimnames = list(names(df.Stack), 
    #                                                          names(df.Stack)))
  cor.matrix <- 1 - abs(stats::cor(df.Stack, method = method))
  dist.matrix <- stats::as.dist(cor.matrix)
  ahc <- stats::hclust(dist.matrix, method = "complete")
  groups <- stats::cutree(ahc, h = 1 - multicollinearity.cutoff)
  if (length(groups) == max(groups)) {
    message(paste("  - No multicollinearity detected in your data at threshold ", 
                  multicollinearity.cutoff, "\n", sep = ""))
    mc <- FALSE
  }
  else {mc <- TRUE }
  if (plot) {
    op <- par(no.readonly = TRUE)
    graphics::par(mar = c(5.1, 5.1, 4.1, 3.1))
    plot(ahc, hang = -1, xlab = "", ylab = "Distance (1 - Pearson's r)", 
         main = "", las = 1, sub = "", axes = F)
    graphics::axis(2, at = seq(0, 1, length = 6), las = 1)
    if (mc) {
      graphics::title(paste("Groups of intercorrelated variables at cutoff", 
                            multicollinearity.cutoff))
      par(xpd = T)
      stats::rect.hclust(ahc, h = 1 - multicollinearity.cutoff)
    }
    else {
      graphics::title(paste("No intercorrelation among variables at cutoff", 
                            multicollinearity.cutoff))
    }
    par(op)
  }
  if (select.variables) {
    sel.vars <- NULL
    for (i in 1:max(groups)) {sel.vars <- c(sel.vars, sample(names(groups[groups == i]), 1))}
  }
  else {
    if (mc) {
      sel.vars <- list()
      for (i in groups) {sel.vars[[i]] <- names(groups)[groups == i]}
    }
    else {sel.vars <- names(df.Stack)}
  }
  return(sel.vars)
}



