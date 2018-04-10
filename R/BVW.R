#' @title Bonus vetus OLS, BVW
#' 
#' @description \code{BVW} estimates gravity models via Bonus 
#' vetus OLS with GDP-weights.
#' 
#' @details \code{Bonus vetus OLS} is an estimation method for gravity models 
#' developed by Baier and Bergstrand (2009, 2010) using GDP-weights to center a
#' Taylor-series (see the references for more information). 
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
#' The \code{BVW} function considers Multilateral Resistance terms and allows to 
#' conduct comparative statics. Country specific effects are subdued due 
#' to demeaning. Hence, unilateral variables apart from \code{inc_o}
#' and \code{inc_d} cannot be included in the estimation.
#' \code{BVW} is designed to be consistent with the Stata code provided at
#' the website
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#' As, to our knowledge at the moment, there is no explicit literature covering 
#' the estimation of a gravity equation by \code{BVW} using panel data, 
#' we do not recommend to apply this method in this case.
#' 
#' @param y name (type: character) of the dependent variable in the dataset 
#' \code{data}, e.g. trade flows. This dependent variable is divided by the 
#' product of unilateral incomes (named \code{inc_o} and \code{inc_d}, e.g. 
#' GDPs or GNPs of the countries of interest) and logged afterwards.
#' The transformed variable is then used as the dependent variable in the 
#' estimation.
#' 
#' @param dist name (type: character) of the distance variable in the dataset 
#' \code{data} containing a measure of distance between all pairs of bilateral
#' partners. It is logged automatically when the function is executed. 
#' 
#' @param x vector of names (type: character) of those bilateral variables in 
#' the dataset \code{data} that should be taken as the independent variables 
#' in the estimation. If an independent variable is a dummy variable,
#' it should be of type numeric (0/1) in the dataset. If an independent variable 
#' is defined as a ratio, it should be logged. Unilateral metric variables 
#' such as GDPs should be inserted via the arguments \code{inc_o} 
#' for the country of origin and \code{inc_d} for the country of destination.
#' As country specific effects are subdued due to demeaning, no further
#' unilateral variables apart from \code{inc_o} and \code{inc_d} can be 
#' added.
#' 
#' @param inc_o variable name (type: character) of the income of the country of 
#' origin in the dataset \code{data}. The dependent variable \code{y} is
#' divided by the product of the incomes \code{inc_d} and \code{inc_o}. 
#' 
#' @param inc_d variable name (type: character) of the income of the country of 
#' destination in the dataset \code{data}. The dependent variable \code{y} is
#' divided by the product of the incomes \code{inc_d} and \code{inc_o}. 
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
#' the estimation of a gravity equation by \code{BVW} 
#' using panel data, cross-sectional data should be used. 
#' 
#' @param ... additional arguments to be passed to \code{BVW}.
#' 
#' @references 
#' For estimating gravity equations via Bonus Vetus OLS see
#' 
#' Baier, S. L. and Bergstrand, J. H. (2009) <DOI:10.1016/j.jinteco.2008.10.004>
#' 
#' Baier, S. L. and Bergstrand, J. H. (2010) in Van Bergeijk, P. A., & Brakman, S. (Eds.) (2010) chapter 4 <DOI:10.1111/j.1467-9396.2011.01000.x>
#' 
#' For more information on gravity models, theoretical foundations and
#' estimation methods in general see 
#' 
#' Anderson, J. E. (1979) <DOI:10.12691/wjssh-2-2-5>
#' 
#' Anderson, J. E. (2010) <DOI:10.3386/w16576>
#' 
#' Anderson, J. E. and van Wincoop, E. (2003) <DOI:10.3386/w8079> 
#' 
#' Head, K., Mayer, T., & Ries, J. (2010) <DOI:10.1016/j.jinteco.2010.01.002>
#' 
#' Head, K. and Mayer, T. (2014) <DOI:10.1016/B978-0-444-54314-1.00003-3>
#' 
#' Santos-Silva, J. M. C. and Tenreyro, S. (2006) <DOI:10.1162/rest.88.4.641> 
#' 
#' and the citations therein.
#' 
#' See \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook} for gravity datasets and Stata code for estimating gravity models.
#' 
#' @examples 
#' \dontrun{
#' data(Gravity_no_zeros)
#' 
#' BVW(y="flow", dist="distw", x=c("rta"), 
#' inc_o="gdp_o", inc_d="gdp_d", vce_robust=TRUE, data=Gravity_no_zeros)
#' 
#' BVW(y="flow", dist="distw", x=c("rta", "comcur", "contig"), 
#' inc_o="gdp_o", inc_d="gdp_d", vce_robust=TRUE, data=Gravity_no_zeros)
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
#' BVW(y="flow", dist="distw", x=c("rta"), inc_o="gdp_o", inc_d="gdp_d", vce_robust=TRUE, data=grav_small)
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
BVW <- function(y, dist, x, inc_o, inc_d, vce_robust=TRUE, data, ...){
  
  if(!is.data.frame(data))                                                stop("'data' must be a 'data.frame'")
  if((vce_robust %in% c(TRUE, FALSE)) == FALSE)                           stop("'vce_robust' has to be either 'TRUE' or 'FALSE'")
  if(!is.character(y)     | !y%in%colnames(data)     | length(y)!=1)      stop("'y' must be a character of length 1 and a colname of 'data'")
  if(!is.character(dist)  | !dist%in%colnames(data)  | length(dist)!=1)   stop("'dist' must be a character of length 1 and a colname of 'data'")
  if(!is.character(x)     | !all(x%in%colnames(data)))                    stop("'x' must be a character vector and all x's have to be colnames of 'data'")  

  if(!is.character(inc_d) | !inc_d%in%colnames(data) | length(inc_d)!=1)  stop("'inc_d' must be a character of length 1 and a colname of 'data'")
  if(!is.character(inc_o) | !inc_o%in%colnames(data) | length(inc_o)!=1)  stop("'inc_o' must be a character of length 1 and a colname of 'data'")
  
  # Transforming data, logging distances ---------------------------------------
  
  d           <- data
  d$dist_log  <- log(d[dist][,1])
  d$count     <- 1:length(d$iso_o)
  
  # Transforming data, logging flows -------------------------------------------
  
  d$y_inc     <- d[y][,1] / (d[inc_o][,1] * d[inc_d][,1])
  d$y_inc_log <- log(d$y_inc)
  
  # GDP weights ----------------------------------------------------------------
  
  inc_world   <- tapply(d[inc_d][,1], d$iso_o, sum)
  d$inc_world <- as.numeric(inc_world[d$iso_o])
  # same for inc_o or inc_d as we have a squared dataset
  
  d$theta_i   <- d[inc_o][,1] / d$inc_world
  d$theta_j   <- d[inc_d][,1] / d$inc_world
  
  # Multilateral resistance (MR) for distance ----------------------------------
  
  d$mr.dist.1 <- NA
  d$mr.dist.2 <- NA
  
  for(i in names(inc_world)){
    
    d[d$iso_o == i,]$mr.dist.1 <- sum(d[d$iso_o == i,]$theta_j * d[d$iso_o == i,]$dist_log)
    d[d$iso_d == i,]$mr.dist.2 <- sum(d[d$iso_d == i,]$theta_i * d[d$iso_d == i,]$dist_log)
  }
  
  d$mr.dist.3   <- sum(d$theta_i * d$theta_j * d$dist_log)
  d$dist_log_mr <- d$dist_log - d$mr.dist.1 - d$mr.dist.2 + d$mr.dist.3
  
  # Multilateral resistance (MR) for the other independent variables -----------
  
  num.ind.var   <- length(x) #independent variables apart from distance
  d_2           <- d
  
  for(j in 1:num.ind.var){
    
    mr.1 <- noquote(paste(c(noquote(x[j]),noquote(".mr1")),collapse="")) 
    mr.2 <- noquote(paste(c(noquote(x[j]),noquote(".mr2")),collapse=""))  
    mr.3 <- noquote(paste(c(noquote(x[j]),noquote(".mr3")),collapse="")) 
    mr   <- noquote(paste(c(noquote(x[j]),noquote(".mr")),collapse="")) 
    
    d_2[mr.1] <- NA
    d_2[mr.2] <- NA
    d_2[mr.3] <- NA
    d_2[mr]   <- NA
    
    for(i in names(inc_world)){
      
      d_2[d_2$iso_o == i,][mr.1] <- sum(d_2[d_2$iso_o == i,]$theta_j * d_2[d_2$iso_o == i,][x[j]])
      d_2[d_2$iso_d == i,][mr.2] <- sum(d_2[d_2$iso_d == i,]$theta_i * d_2[d_2$iso_d == i,][x[j]])
    }
    
    d_2[mr.3] <- sum(d_2$theta_i * d_2$theta_j * d_2[x[j]])
    d_2[x[j]] <- d_2[x[j]] - d_2[mr.1] - d_2[mr.2] + d_2[mr.3]
    
  }
  
  # Model ----------------------------------------------------------------------
  
  x_mr <- paste0(x,"_mr")
  
  # new row in dataset for independent _mr variable
  for(j in x){
    l       <- which(x == j)
    mr      <- x_mr[l]
    d_2[mr] <- NA
    d_2[mr] <- d_2[x[l]]
  }
  
  vars      <- paste(c("dist_log_mr", x_mr), collapse=" + ")
  form      <- paste("y_inc_log","~",vars)
  form2     <- stats::as.formula(form)
  model.BVW <- stats::lm(form2, data = d_2)
  
  # Return ---------------------------------------------------------------------
  
  if(vce_robust == TRUE){
    return.object.1      <- .robustsummary.lm(model.BVW, robust=TRUE)
    return.object.1$call <- form2
    return(return.object.1)}
  
  if(vce_robust == FALSE){
    return.object.1      <- .robustsummary.lm(model.BVW, robust=FALSE)
    return.object.1$call <- form2
    return(return.object.1)}
}