#' @title Double Demeaning, DDM
#' 
#' @description \code{DDM} estimates gravity models via double demeaning the 
#' left hand side and right hand side of the gravity equation.
#' 
#' @details \code{DDM} is an estimation method for gravity models presented
#' in Head and Mayer (2014) (see the references for more information).  
#' To execute the function a square gravity dataset with all pairs of 
#' countries, ISO-codes for the country of origin and destination, a measure of 
#' distance between the bilateral partners as well as all 
#' information that should be considered as dependent an independent 
#' variables is needed. 
#' Make sure the ISO-codes are of type "character".
#' Missing bilateral flows as well as incomplete rows should be 
#' excluded from the dataset. 
#' Furthermore, flows equal to zero should be excluded as the gravity equation 
#' is estimated in its additive form.  
#' Country specific effects are subdued due double demeaning. 
#' Hence, unilateral income proxies such as GDP cannot be 
#' considered as exogenous variables.
#' 
#' \code{DDM} is designed to be consistent with the Stata code provided at
#' the website
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#' As, to our knowledge at the moment, there is no explicit literature covering 
#' the estimation of a gravity equation by \code{DDM} using panel data, 
#' we do not recommend to apply this method in this case.
#' 
#' @param y name (type: character) of the dependent variable in the dataset 
#' \code{data}, e.g. trade flows. This variable is logged and then used as the 
#' dependent variable in the estimation.
#' 
#' @param dist name (type: character) of the distance variable in the dataset 
#' \code{data} containing a measure of distance between all pairs of bilateral
#' partners. It is logged automatically when the function is executed. 
#' 
#' @param x vector of names (type: character) of those bilateral variables in 
#' the dataset \code{data} that should be taken as the independent variables 
#' in the estimation. If an independent variable is a dummy variable,
#' it should be of type numeric (0/1) in the dataset. If an independent variable 
#' is defined as a ratio, it should be logged. 
#' Unilateral effects drop out due to double demeaning and therefore 
#' cannot be estimated.
#' 
#' @param vce_robust robust (type: logic) determines whether a robust 
#' variance-covariance matrix should be used. The default is set to \code{TRUE}. 
#' If set \code{TRUE} the estimation results are consistent with the 
#' Stata code provided at the website
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#' 
#' @param data name of the dataset to be used (type: character). 
#' To estimate gravity equations, a square gravity dataset including bilateral 
#' flows defined by the argument \code{y}, ISO-codes of type character 
#' (called \code{iso_o} for the country of origin and \code{iso_d} for the 
#' destination country), a distance measure defined by the argument \code{dist} 
#' and other potential influences given as a vector in \code{x} are required. 
#' All dummy variables should be of type numeric (0/1). Missing trade flows as 
#' well as incomplete rows should be excluded from the dataset. 
#' Furthermore, flows equal to zero should be excluded as the gravity equation 
#' is estimated in its additive form.
#' As, to our knowledge at the moment, there is no explicit literature covering 
#' the estimation of a gravity equation by \code{DDM} 
#' using panel data, cross-sectional data should be used. 
#' 
#' @param ... additional arguments to be passed to \code{DDM}.
#' 
#' @references 
#' For more information on Double Demeaning as well as information on gravity 
#' models, theoretical foundations and estimation methods in general see
#' 
#' Head, K. and Mayer, T. (2014) <DOI:10.1016/B978-0-444-54314-1.00003-3>
#' 
#' as well as
#' 
#' Anderson, J. E. (1979) <DOI:10.12691/wjssh-2-2-5>
#' 
#' Anderson, J. E. (2010) <DOI:10.3386/w16576>
#' 
#' Anderson, J. E. and van Wincoop, E. (2003) <DOI:10.3386/w8079> 
#' 
#' Baier, S. L. and Bergstrand, J. H. (2009) <DOI:10.1016/j.jinteco.2008.10.004>
#' 
#' Baier, S. L. and Bergstrand, J. H. (2010) in Van Bergeijk, P. A., & Brakman, S. (Eds.) (2010) chapter 4 <DOI:10.1111/j.1467-9396.2011.01000.x>
#' 
#' Head, K., Mayer, T., & Ries, J. (2010) <DOI:10.1016/j.jinteco.2010.01.002>
#' 
#' Santos-Silva, J. M. C. and Tenreyro, S. (2006) <DOI:10.1162/rest.88.4.641> 
#' 
#' and the citations therein.
#' 
#' 
#' See \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook} for gravity datasets and Stata code for estimating gravity models.
#' 
#' @examples 
#' \dontrun{
#' data(Gravity_no_zeros)
#' 
#' DDM(y="flow", dist="distw", x=c("rta"), 
#' vce_robust=TRUE, data=Gravity_no_zeros)
#' 
#' DDM(y="flow", dist="distw", x=c("rta", "comcur", "contig"), 
#' vce_robust=TRUE, data=Gravity_no_zeros)
#' }
#' 
#' \dontshow{
#' # examples for CRAN checks:
#' # executable in < 5 sec together with the examples above
#' # not shown to users
#' 
#' data(Gravity_no_zeros)
#' # choose exemplarily 10 biggest countries for check data
#' countries_chosen <- names(sort(table(Gravity_no_zeros$iso_o), decreasing = TRUE)[1:10])
#' grav_small <- Gravity_no_zeros[Gravity_no_zeros$iso_o %in% countries_chosen,]
#' DDM(y="flow", dist="distw", x=c("rta"), vce_robust=TRUE, data=grav_small)
#' }
#' 
#' @return
#' The function returns the summary of the estimated gravity model as an 
#' \code{\link[stats]{lm}}-object.
#' 
#' @seealso \code{\link[stats]{lm}}, \code{\link[lmtest]{coeftest}}, 
#' \code{\link[sandwich]{vcovHC}}
#' 
#' @export 
#' 
DDM <- function(y, dist, x, vce_robust=TRUE, data, ...){
  
  if(!is.data.frame(data))                                                stop("'data' must be a 'data.frame'")
  if((vce_robust %in% c(TRUE, FALSE)) == FALSE)                           stop("'vce_robust' has to be either 'TRUE' or 'FALSE'")
  if(!is.character(y)     | !y%in%colnames(data)     | length(y)!=1)      stop("'y' must be a character of length 1 and a colname of 'data'")
  if(!is.character(dist)  | !dist%in%colnames(data)  | length(dist)!=1)   stop("'dist' must be a character of length 1 and a colname of 'data'")
  if(!is.character(x)     | !all(x%in%colnames(data)))                    stop("'x' must be a character vector and all x's have to be colnames of 'data'")  

  # Transforming data, logging distances ---------------------------------------
  
  d               <- data
  d$dist_log      <- (log(d[dist][,1]))
  d$count         <- 1:length(d$iso_o)
  
  # Transforming data, logging flows -------------------------------------------
  
  d$y_log         <- log(d[y][,1])
  
  # Substracting the means -----------------------------------------------------
  
  d$y_log_dd      <- rep(NA, times = length(d$dist_log))
  d$dist_log_dd   <- rep(NA, times = length(d$dist_log))
  
  mean.y_log.1    <- tapply(d$y_log, d$iso_o, mean)
  mean.y_log.2    <- tapply(d$y_log, d$iso_d, mean)  
  
  mean.dist_log.1 <- tapply(d$dist_log, d$iso_o, mean)
  mean.dist_log.2 <- tapply(d$dist_log, d$iso_d, mean)

  d$y_log_dd      <- d$y_log
  d$dist_log_dd   <- d$dist_log
  
  for(i in unique(d$iso_o)){
    d[d$iso_o == i,]$y_log_dd    <- d[d$iso_o == i,]$y_log_dd - mean.y_log.1[i]
    d[d$iso_o == i,]$dist_log_dd <- d[d$iso_o == i,]$dist_log_dd - mean.dist_log.1[i]
  }
  
  for(i in unique(d$iso_d)){
    d[d$iso_d == i,]$y_log_dd    <- d[d$iso_d == i,]$y_log_dd - mean.y_log.2[i]
    d[d$iso_d == i,]$dist_log_dd <- d[d$iso_d == i,]$dist_log_dd - mean.dist_log.2[i]
  }
  
  d$y_log_dd    <- d$y_log_dd + mean(d$y_log)
  d$dist_log_dd <-  d$dist_log_dd + mean(d$dist_log)
  
  # Substracting the means for the other independent variables -----------------
  
  ind.var.dd     <- list(length=length(x))
  mean.ind.var.1 <- list(legth=length(x))
  mean.ind.var.2 <- list(legth=length(x))
  
  for(j in 1:length(x)){
    ind.var.dd[[j]]     <- d[x[j]][,1]
    mean.ind.var.1[[j]] <- tapply(d[x[j]][,1], d$iso_o, mean)
    mean.ind.var.2[[j]] <- tapply(d[x[j]][,1], d$iso_d, mean)
  }
  
  w   <- letters[1:length(x)]
  
  d_2 <- d
  for(j in 1:length(x)){
    d_2[w[j]] <- ind.var.dd[[j]]
  }

  d_3 <- d_2
  
  for(j in 1:length(x)){
    
    for(i in unique(d_2$iso_o)){
      d_2[d_2$iso_o == i,][w[j]] <- d_2[d_2$iso_o == i,][w[j]] - mean.ind.var.1[[j]][i]
    }
    
    for(i in unique(d_2$iso_d)){
      d_2[d_2$iso_d == i,][w[j]] <- d_2[d_2$iso_d == i,][w[j]] - mean.ind.var.2[[j]][i]
    }
    
    d_2[w[j]] <- d_2[w[j]] + mean(d_2[x[j]][,1])
    d_3[x[j]] <- d_2[w[j]]
  }
  
  # Model ----------------------------------------------------------------------
  
  x_dd <- paste0(x,"_dd")
  
  # new row in dataset for independent _dd variable
  for(j in x){
    l       <- which(x == j)
    dd      <- x_dd[l]
    d_3[dd] <- NA
    d_3[dd] <- d_3[x[l]]
  }
  
  vars      <- paste(c("dist_log_dd", x_dd), collapse=" + ")
  form      <- paste("y_log_dd", "~", vars, "- 1")
  form2     <- stats::as.formula(form)
  
  model.DDM <- stats::lm(form2, data = d_3)   

  # Return ---------------------------------------------------------------------
  
  if(vce_robust == TRUE){
    return.object.1         <- .robustsummary.lm(model.DDM, robust=TRUE)
    return.object.1$call    <- form2
    return(return.object.1)}
  
  if(vce_robust == FALSE){
    return.object.1        <- .robustsummary.lm(model.DDM, robust=FALSE)
    return.object.1$call   <- form2
    return(return.object.1)}
}