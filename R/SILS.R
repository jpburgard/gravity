#' @title Structural Iterated Least Squares, SILS
#' 
#' @description \code{SILS} estimates gravity models via 
#' Structural Iterated Least Squares and an explicit inclusion
#' of the Multilateral Resistance terms.
#' 
#' @details \code{SILS} is an estimation method for gravity models
#' developed by Head and Mayer (2014) (see the references for
#' more information). 
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
#' The function \code{SILS} utilizes the relationship between the Multilateral 
#' Resistance terms and the transaction costs. The parameters are estimated by 
#' an iterative procedure. The function executes loops until the parameters 
#' stop changing significantly.
#' 
#' \code{SILS} is designed to be consistent with the Stata code provided at
#' the website
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#' As, to our knowledge at the moment, there is no explicit literature covering 
#' the estimation of a gravity equation by \code{SILS} using panel data, 
#' we do not recommend to apply this method in this case.
#' 
#' @param y name (type: character) of the dependent variable in the dataset 
#' \code{data}, e.g. trade flows. This dependent variable is divided by the 
#' product of unilateral incomes (\code{inc_o} and \code{inc_d}, e.g. 
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
#' 
#' @param inc_o variable name (type: character) of the income of the country of 
#' origin in the dataset \code{data}. The dependent variable \code{y} is
#' divided by the product of the incomes \code{inc_d} and \code{inc_o}. 
#' 
#' @param inc_d variable name (type: character) of the income of the country of 
#' destination in the dataset \code{data}. The dependent variable \code{y} is
#' divided by the product of the incomes \code{inc_d} and \code{inc_o}. 
#' 
#' @param maxloop maximum number of outer loop iterations. 
#' The default is set to 50. There will be a warning if the iterations
#' did not converge.
#' 
#' @param maxloop2 maximum number of inner loop iterations. 
#' The default is set to 50. There will be a warning if the iterations
#' did not converge.
#' 
#' @param dec_places number of decimal places that should not change after a new 
#' iteration for the estimation to stop. The default is set to 4. 
#' 
#' @param vce_robust robust (type: logic) determines whether a robust 
#' variance-covariance matrix should be used. The default is set to \code{TRUE}. 
#' If set \code{TRUE} the estimation results are consistent with the 
#' Stata code provided at the website
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#' 
#' @param verbose (type: logic) determines whether the estimated coefficients
#' of each iteration should be printed in the console. The default is set
#' to \code{FALSE}.
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
#' the estimation of a gravity equation by \code{SILS} 
#' using panel data, cross-sectional data should be used. 
#' 
#' @param ... additional arguments to be passed to functions used by 
#' \code{SILS}.
#' 
#' @references 
#' For information on \code{SILS} as well as more information on gravity 
#' models, theoretical foundations and suitable estimation methods in general see
#' 
#' Head, K. and Mayer, T. (2014) <DOI:10.1016/B978-0-444-54314-1.00003-3>
#' 
#' and
#' 
#' Anderson, J. E. and van Wincoop, E. (2003) <DOI:10.3386/w8079> 
#' 
#' as well as
#' 
#' Anderson, J. E. (1979) <DOI:10.12691/wjssh-2-2-5>
#' 
#' Anderson, J. E. (2010) <DOI:10.3386/w16576>
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
#' SILS(y="flow", dist="distw", x=c("rta"), inc_o="gdp_o", inc_d="gdp_d", 
#' maxloop=50, maxloop2=50, dec_places=4, vce_robust=TRUE, verbose=FALSE, 
#' data=Gravity_no_zeros)
#' 
#' SILS(y="flow", dist="distw", x=c("rta", "comcur", "contig"), 
#' inc_o="gdp_o", inc_d="gdp_d", maxloop=50, maxloop2=50, dec_places=4, 
#' vce_robust=TRUE, verbose=TRUE, data=Gravity_no_zeros)
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
#' SILS(y="flow", dist="distw", x=c("rta"), inc_o="gdp_o", inc_d="gdp_d", maxloop=50, maxloop2=50, dec_places=4, vce_robust=TRUE, verbose=TRUE, data=grav_small)
#' }
#' 
#' @return
#' The function returns the summary of the estimated gravity model as an 
#' \code{\link[stats]{lm}}-object. It furthermore returns the resulting coefficients for each
#' iteration.
#' 
#' @seealso \code{\link[stats]{lm}}, \code{\link[lmtest]{coeftest}}, 
#' \code{\link[sandwich]{vcovHC}}
#' 
#' 
#' @export 
#' 
SILS <- function(y, dist, x, inc_o, inc_d, maxloop=50, maxloop2=50, dec_places=4, vce_robust=TRUE, verbose=FALSE, data, ...){
  
  if(!is.data.frame(data))                                                stop("'data' must be a 'data.frame'")
  if((vce_robust %in% c(TRUE, FALSE)) == FALSE)                           stop("'vce_robust' has to be either 'TRUE' or 'FALSE'")
  if(!is.character(y)     | !y%in%colnames(data)     | length(y)!=1)      stop("'y' must be a character of length 1 and a colname of 'data'")
  if(!is.character(dist)  | !dist%in%colnames(data)  | length(dist)!=1)   stop("'dist' must be a character of length 1 and a colname of 'data'")
  if(!is.character(x)     | !all(x%in%colnames(data)))                    stop("'x' must be a character vector and all x's have to be colnames of 'data'")  

  if(!is.character(inc_d) | !inc_d%in%colnames(data) | length(inc_d)!=1)  stop("'inc_d' must be a character of length 1 and a colname of 'data'")
  if(!is.character(inc_o) | !inc_o%in%colnames(data) | length(inc_o)!=1)  stop("'inc_o' must be a character of length 1 and a colname of 'data'")
  if(maxloop!=round(maxloop) | maxloop2!=round(maxloop2) | maxloop<=0 | maxloop2<=0) stop("'maxloop' and 'maxloop2' have to be integers with values >0")
  
  # checking whether parameters are valid --------------------------------------
  
  if(maxloop < 1)                                                         stop("maxloop has to be an integer greater or equal 1")
  if(maxloop2 < 1)                                                        stop("maxloop2 has to be an integer greater or equal 1")
  if(dec_places < 1)                                                      stop("dec_places should equal 1 or a higher number")
  
  d                <- data
  d$dist_log       <- (log(d[dist][,1]))
  
  # Setting starting values for the first iteration ----------------------------
  
  d$P_i            <- 1
  d$P_j            <- 1
  
  loop             <- 0
  dec.point        <- 1*10^-dec_places
  beta_dist        <- 1
  beta_dist_old    <- 0
  coef.dist        <- 1
  
  beta             <- vector(length=length(x))
  names(beta)      <- x
  for(j in 1:length(x)){beta[j] <- 1}
  beta_old         <- vector(length=length(x))
  names(beta_old)  <- x
  
  for(j in 1:length(x)){beta_old[j] <- 0}
  
  coef_x           <- data.frame(matrix(nrow=1, ncol=length(x)))
  coef_x[1,]       <- 1
  names(coef_x)    <- x
  
  # Begin iterations -----------------------------------------------------------
  
  while(loop <= maxloop & 
        abs(beta_dist - beta_dist_old) > dec.point & 
        prod(abs(beta - beta_old) > dec.point) == 1){
    
    # Updating betas -----------------------------------------------------------
    
    beta_dist_old <- beta_dist
    for(j in 1:length(x)){beta_old[j] <- beta[j]}
    
    # Updating transaction costs -----------------------------------------------
    
    costs_a <- data.frame(matrix(nrow=length(d[y][,1]), ncol=length(x)))
    for(j in 1:length(x)){costs_a[,j] <- beta[j] * d[x[j]][,1]}
    costs_b <- apply(X=costs_a, MARGIN=1, FUN=sum)
    d$t_ij  <- exp(beta_dist * d$dist_log + costs_b) 
    
    # Contraction mapping ------------------------------------------------------
    
    d$P_j_old <- 0
    d$P_i_old <- 0
    j         <- 1
    
    while(j <= maxloop2 & 
          sum(abs(d$P_j - d$P_j_old)) > dec.point & 
          sum(abs(d$P_i - d$P_i_old)) > dec.point){
      
      d$P_j_old <- d$P_j
      d$P_i_old <- d$P_i
      
      # Inward MR --------------------------------------------------------------
      
      for(i in unique(d$iso_d)){
        d$P_j[d$iso_d==i] <- sum(d$t_ij[d$iso_d==i] * 
                                   d[inc_o][,1][d$iso_d==i] / 
                                   d$P_i[d$iso_d==i])
      }
      
      # Outwad MR --------------------------------------------------------------
      
      for(i in unique(d$iso_o)){
        d$P_i[d$iso_o==i] <- sum(d$t_ij[d$iso_o==i] * 
                                   d[inc_d][,1][d$iso_o==i] / 
                                   d$P_j[d$iso_o==i])}
      
      j <- j+1
      
      if(j == maxloop2){
        warning("The inner iteration did not converge before the inner loop reached maxloop2=",maxloop2," iterations")
      }
      
    }
    
    # Model --------------------------------------------------------------------
    
    vars       <- paste(c("dist_log", x), collapse=" + ")
    form       <- paste("log(d[y][,1]) - log((d[inc_o][,1] * d[inc_d][,1]) /
                    (d$P_i * d$P_j))","~",vars)
    form2      <- stats::as.formula(form)
    
    model.SILS <- stats::lm(form2, data = d)
    
    # Updating coefficients ----------------------------------------------------
    
    beta_dist <- stats::coef(model.SILS)[2]
    for(j in 1:length(x)){beta[j] <- stats::coef(model.SILS)[j+2]}
    
    coef.dist <- c(coef.dist, beta_dist)
    coef_x    <- rbind(coef_x, rep(0, times=length(x)))
    for(j in 1:length(x)){coef_x[x[j]][loop+2,] <- beta[j]}
    
    # Coefficients -------------------------------------------------------------
    
    coef.SILS <- cbind(loop=c(1:loop), dist=as.numeric(coef.dist)[2:(loop+1)], 
                       coef_x[2:(loop+1),]) 
    
    if(verbose == TRUE){
      cat("This is round",loop,"\n")
      cat("The coefficients are",beta_dist,beta,"\n")
    }
    
    loop <- loop + 1 
    
    if(loop == maxloop){
      warning("The outer iteration did not converge before the outer loop reached maxloop=",maxloop," iterations")
    }
    
  }
  
  # Return --------------------------------------------------------------------- 
  
  if(vce_robust == TRUE){
    return.object.1      <- .robustsummary.lm(model.SILS, robust=TRUE)
    return.object.1$call <- form2
    return(return.object.1)}
  
  if(vce_robust == FALSE){
    return.object.1      <- .robustsummary.lm(model.SILS, robust=FALSE)
    return.object.1$call <- form2
    return(return.object.1)}
  
}