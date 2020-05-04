#' Variable Selection
#'
#' Returns plot with correlation values among predicted variables.
#'
varSelection <- function(df, 
                         VarMethod = NULL, # c('vif','correlation','both'), 
                         vifs=NULL, 
                         cutoff=NULL,
                         removeVariable=NULL,
                         keep_var=FALSE, ...)
  
#'  Variable selection using two approaches: vif and person correlation.
#'
#' @param df dataframe-class. 
#' @param VarMethod Method used to variable selection. \code{vif} uses variance inflation factor 
#' based on \code{car} package. \code{correlation} uses pearson correlation to 
#' select random variable using cutoff parameter. \code{both} use a two approaches describe above.
#' @param vifs numeric vector. vif uses numeric vector to return variables selection.
#' @param cutoff numeric vector. Person correlation coef. used like threshold.
#' @param removeVariable Variables choosed to remove before running varSelection.
#' @param keep_var logical. If \code{TRUE}, all variables will be save in \code{Variables_Selected} folder.
#' (\code{default = FALSE}).
#' @param ... no implemented. Use \code{trace=TRUE} to show vif running. \code{default=FALSE}.
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
{
  
  message('Variable selection...')
  
  if(!class(df)=="data.frame") {
    stop(paste0('df must be a dataframe object'))
  }
    
  if(is.null(VarMethod)){
    stop(paste0('Choose VarMethod between vif, correlation, or both'))
  } 
  
  if(VarMethod=='vif' & is.null(vifs)){
    stop(paste0('vifs parameter must be a number or numeric vector e.g. vifs=2 or vifs=c(2,5)'))
  } 
  if(VarMethod=='correlation'){
    if(is.null(cutoff)) {stop(paste0('cutoff parameter must be a number or numeric vector e.g. cutoff=2 or cutoff=c(2,5)'))}
    if(cutoff <= 0 | cutoff >=1 ){stop(paste('cutoff parameter must be a number from 0 to 1'))}
  }
  
  if(VarMethod=='both'){
   if(is.null(vifs) | is.null(cutoff)){
     stop(paste0('You need to choose a vif AND cutoff parameter'))
   } 
  }
  
 
  ValsENMV <- df
  
  if(!is.null(removeVariable)){
    Dataset <-  as.data.frame(ValsENMV) %>% dplyr::select(-one_of(removeVariable)) %>%
      as.matrix()
  } else{Dataset <- ValsENMV}
  
  if(!class(Dataset) =='data.frame') {message('DataSet must be data.frame object')}
  
  if(VarMethod=='vif'){
    if(is.null(vifs)){message('We a number or a vector of numbers')}
    vifNames <- paste0('vif', vifs)
    vars.Select <- list()
    for(i in 1:length(vifs)){
      VIFs <- VIF_Selection(df.Stack=Dataset, thresh=vifs[i],...)
      vars.Select[[i]] <- VIFs
    }
    namesSets <- vifNames
  }
  
  if(VarMethod=='correlation'){
    corNames <- paste0('cor', cutoff)
    vars.Select <- list()
    for(i in 1:length(cutoff)){
      vars.Select[[i]] <- MTHs_Heuristico(df.Stack=Dataset, multicollinearity.cutoff = cutoff[i]) 
    }
    
    namesSets <- corNames
  }
  
  
  if(VarMethod=='both'){
    vifNames <- paste0('vif', vifs)
    corNames <- paste0('cor', cutoff)
    
    vars.m1 <- list()
    for(i in 1:length(vifs)){
      vars.m1[[i]] <- VIF_Selection(df.Stack=Dataset, thresh=vifs[i])
    }
    vars.m2<- list()
    for(i in 1:length(cutoff)){
      vars.m2[[i]] <- MTHs_Heuristico(df.Stack=Dataset, multicollinearity.cutoff = cutoff[i]) 
    }
    vars.Select <- do.call(c, list(vars.m1, vars.m2))
    namesSets <- c(vifNames, corNames)
  }
  
  VarSel <- new("VarSelection", variables = vars.Select, VarsSets = namesSets)
  
  if(!keep_var==FALSE){
    dir.create('Variables_Selected')
    for(f in 1: length(VarSel@variables)){
      varSetsExport <- data.frame(variables=VarSel@variables[[f]], CODE=VarSel@VarsSets[f])
      
      
      readr::write_csv(varSetsExport, paste0('./Variables_Selected/Var_Set_',f,'.csv'))
    } 
  }
  VarSel
  message(paste0(length(VarSel@variables), ' variable selected'))
  return(VarSel)
  
  
  
}
