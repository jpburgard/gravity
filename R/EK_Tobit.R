#' @title Eaton and Kortum (2001) Tobit model, EK_Tobit
#' 
#' @description \code{EK_Tobit} estimates gravity models in their additive form
#' by conducting a censored regression.
#' It follows the Eaton and Kortum (2001) Tobit model where each country 
#' is assigned specific ceonsoring bounds.
#' 
#' @details \code{EK_Tobit} represents the Eaton and Kortum (2001) Tobit model.
#' When taking the log of the gravity equation flows equal to zero
#' constitute a problem as their log is not defined.
#' Therefore, in \code{EK_Tobit} all values of the dependent variable 
#' are redefined as intervals.
#' The positive observations have both interval bounds equal
#' to their original value. 
#' For zero flows the interval is left open. The right border
#' of the interval is set to the log of the minimum positive trade flow of 
#' the respective importing country.
#' The defined data object of class \code{\link[survival]{Surv}} is
#' then inserted in \code{\link[survival]{survreg}} for the 
#' parameter estimation.
#'  
#' To execute the function a square gravity dataset with all pairs of 
#' countries, ISO-codes for the country of origin and destination, a measure of 
#' distance between the bilateral partners as well as all 
#' information that should be considered as dependent an independent 
#' variables is needed. 
#' Missing bilateral flows as well as incomplete rows should be 
#' excluded from the dataset.  
#' Zero trade flows are allowed. 
#' 
#' \code{EK_Tobit} is designed to be consistent with the Stata code provided at
#' the website
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#' 
#' Up to now, the function is designed for cross-sectional data,
#' but can be extended to panel data using the 
#' \code{\link[survival]{survreg}} function.
#' 
#' For other Tobit functions, see \code{\link[gravity]{Tobit}}
#' for a simple Tobit model where number \code{1} is added to all observations
#' and \code{\link[gravity]{ET_Tobit}} for the Eaton and Tamura (1994) 
#' threshold Tobit model where instead of simply adding number \code{1} 
#' to the data the threshold is estimated.
#' 
#' @param y name (type: character) of the dependent variable in the dataset 
#' \code{data}, e.g. trade flows.
#' The variable is logged and then taken as the dependent variable in 
#' the regression. As the log of zero is not defined, 
#' all flows equal to zero are replaced by
#' a left open interval with the logged minimum trade flow of the
#' respective importing country as right border.
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
#' Unilateral variables such as country dummies or incomes can be added. 
#' If unilateral metric variables such as GDPs should be used as independent 
#' variables, those variables have to be logged first and the 
#' logged variable can be used in \code{x}.
#' Interaction terms can be added.
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
#' Zero trade flows are allowed.
#' 
#' @param ... additional arguments to be passed to \code{EK_Tobit}.
#' 
#' @references 
#' For more information on gravity models, theoretical foundations and
#' estimation methods in general see 
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
#' Head, K. and Mayer, T. (2014) <DOI:10.1016/B978-0-444-54314-1.00003-3>
#' 
#' Santos-Silva, J. M. C. and Tenreyro, S. (2006) <DOI:10.1162/rest.88.4.641> 
#' 
#' and the citations therein.
#' 
#' 
#' Especially for Tobit models see
#' 
#' Tobin, J. (1958) <DOI:10.2307/1907382>
#' 
#' Eaton, J., & Tamura, A. (1994) <DOI:10.3386/w4758>
#' 
#' Eaton, J., & Kortum, S. (2001) <DOI:10.3386/w8070>.
#' 
#' 
#' See Carson, R. T., & Sun, Yixiao (2007) <DOI:10.1111/j.1368-423X.2007.00218.x>
#' for the estimation of the threshold in a Tobit model. 
#' 
#' See \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook} for gravity datasets and Stata code for estimating gravity models.
#' 
#' @examples 
#' \dontrun{
#' # Example for data with zero trade flows
#' data(Gravity_zeros)
#' 
#' Gravity_zeros$lgdp_o <- log(Gravity_zeros$gdp_o)
#' Gravity_zeros$lgdp_d <- log(Gravity_zeros$gdp_d)
#' 
#' EK_Tobit(y="flow", dist="distw", x=c("rta","lgdp_o","lgdp_d"), 
#' vce_robust=TRUE, data=Gravity_zeros)
#' 
#' EK_Tobit(y="flow", dist="distw", x=c("rta","iso_o","iso_d"), 
#' vce_robust=TRUE, data=Gravity_zeros)
#' }
#' 
#' \dontshow{
#' # examples for CRAN checks:
#' # executable in < 5 sec together with the examples above
#' # not shown to users
#' 
#' data(Gravity_zeros)
#' Gravity_zeros$lgdp_o <- log(Gravity_zeros$gdp_o)
#' Gravity_zeros$lgdp_d <- log(Gravity_zeros$gdp_d)
#' 
#' # choose exemplarily 10 biggest countries for check data
#' countries_chosen_zeros <- names(sort(table(Gravity_zeros$iso_o), decreasing = TRUE)[1:10])
#' grav_small_zeros <- Gravity_zeros[Gravity_zeros$iso_o %in% countries_chosen_zeros,]
#' EK_Tobit(y="flow", dist="distw", x=c("rta","lgdp_o","lgdp_d"), vce_robust=TRUE, data=grav_small_zeros)
#' }
#' 
#' @return
#' The function returns the summary of the estimated gravity model as a 
#' \code{\link[survival]{survreg}}-object.
#' 
#' @seealso \code{\link[survival]{Surv}}, \code{\link[survival]{survreg}}
#' 
#' @export 
#' 
EK_Tobit <- function(y, dist, x, vce_robust=TRUE, data, ...){
  
  if(!is.data.frame(data))                                                stop("'data' must be a 'data.frame'")
  if((vce_robust %in% c(TRUE, FALSE)) == FALSE)                           stop("'vce_robust' has to be either 'TRUE' or 'FALSE'")
  if(!is.character(y)     | !y%in%colnames(data)     | length(y)!=1)      stop("'y' must be a character of length 1 and a colname of 'data'")
  if(!is.character(dist)  | !dist%in%colnames(data)  | length(dist)!=1)   stop("'dist' must be a character of length 1 and a colname of 'data'")
  if(!is.character(x)     | !all(x%in%colnames(data)))                    stop("'x' must be a character vector and all x's have to be colnames of 'data'")  
  
  # Transforming data, logging flows, distances --------------------------------

  d                           <- data
  d$dist_log                  <- (log(d[dist][,1]))
  d$y_log                     <- log(d[y][,1]) 
  d[d[y][,1] == 0,]$y_log     <- NA
  
  # Setting the lower and upper bounds -----------------------------------------
  
  # lower bound
  d$l_flows_EK1               <- NA
  
  # upper bound
  d$l_flows_EK2               <- NA
  
  # if trade flows > 0: lower bound = upper bound = logged flows
  d[d$flow > 0,]$l_flows_EK1  <- d[d[y][,1] > 0,]$l_flows_EK2 <- d[d[y][,1] > 0,]$y_log
  
  # if trade flows == 0
  d[d$flow == 0,]$l_flows_EK1 <- -Inf 
  
  # we need ido_o and iso_d
  exportmin                   <- sapply(unique(d["iso_d"][,1]), 
                                        function(iso) 
                                          min(d[d["iso_d"][,1] == iso,]$y_log, 
                                              na.rm = TRUE))
  
  for(iso in unique(d["iso_d"][,1])){
    # iso = unique(d["iso_d"][,1])[1]
    if(nrow(d[d$flow == 0 & d["iso_d"][,1] == iso,]) > 0){
      d[d$flow == 0 & d["iso_d"][,1] == iso,]$l_flows_EK2 <- exportmin[iso]
    }
  }
  
  # Create Survival Object
  dta_int              <- survival::Surv(d$l_flows_EK1, d$l_flows_EK2, type = "interval2")
  
  # Model ----------------------------------------------------------------------
  
  vars                 <- paste(c("dist_log", x), collapse = " + ")
  form                 <- paste("dta_int", "~", vars)
  form2                <- stats::as.formula(form)
  model.EK_Tobit       <- survival::survreg(form2, data = d, dist = "gaussian", robust = vce_robust)
  
  # Return --------------------------------------------------------------------- 
  
  return.object.1      <- summary(model.EK_Tobit)
  return.object.1$call <- form2
  return(return.object.1)
}

