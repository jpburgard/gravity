#' @title Negative Binomial Pseudo Maximum Likelihood, NBPML
#' 
#' @description \code{NBPML} estimates gravity models in their 
#' multiplicative form via Negative Binomial Pseudo Maximum Likelihood.
#' 
#' @details \code{NBPML} is an estimation method for gravity models
#' belonging to generalized linear models.
#' It is estimated via \code{\link[MASS]{glm.nb}}  
#' using the negative binomial distribution and a log-link.
#' To execute the function a square gravity dataset with all pairs of 
#' countries, ISO-codes for the country of origin and destination, a measure of 
#' distance between the bilateral partners as well as all 
#' information that should be considered as dependent an independent 
#' variables is needed. 
#' Make sure the ISO-codes are of type "character".
#' Missing bilateral flows as well as incomplete rows should be 
#' excluded from the dataset.  
#' Zero trade flows are allowed.
#' For similar functions, utilizing the multiplicative form via the log-link, 
#' but different distributions, see \code{\link[gravity]{PPML}}, \code{\link[gravity]{GPML}}, and \code{\link[gravity]{NLS}}.
#' 
#' \code{NBPML} estimation can be used for both, cross-sectional as well as 
#' panel data. 
#' It is up to the user to ensure that the functions can be applied 
#' to panel data. 
#' Depending on the panel dataset and the variables - 
#' specifically the type of fixed effects - 
#' included in the model, it may easily occur that the model is not computable. 
#' Also, note that by including bilateral fixed effects such as country-pair 
#' effects, the coefficients of time-invariant observables such as distance 
#' can no longer be estimated. 
#' Depending on the specific model, the code of the 
#' respective function may has to be changed in order to exclude the distance 
#' variable from the estimation. 
#' At the very least, the user should take special 
#' care with respect to the meaning of the estimated coefficients and variances 
#' as well as the decision about which effects to include in the estimation. 
#' When using panel data, the parameter and variance estimation of the models 
#' may have to be changed accordingly.
#' For a comprehensive overview of gravity models for panel data 
#' see Egger and Pfaffermayr (2003), Gomez-Herrera (2013) and Head, Mayer and 
#' Ries (2010) as well as the references therein. 
#' 
#' @param y name (type: character) of the dependent variable in the dataset 
#' \code{data}, e.g. trade flows. 
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
#' When using panel data, a variable for the time may be included in the 
#' dataset. Note that the variable for the time dimension should be of 
#' type: factor. See the references for more information on panel data.
#' 
#' @param ... additional arguments to be passed to functions used by 
#' \code{NBPML}.
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
#' See \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook} for gravity datasets and Stata code for estimating gravity models.
#' 
#' 
#' For estimating gravity equations using panel data see 
#' 
#' Egger, P., & Pfaffermayr, M. (2003) <DOI:10.1007/s001810200146>
#' 
#' Gomez-Herrera, E. (2013) <DOI:10.1007/s00181-012-0576-2>
#' 
#' and the references therein.
#' 
#' @examples 
#' \dontrun{
#' # Example for data with zero trade flows
#' data(Gravity_zeros)
#' 
#' NBPML(y="flow", dist="distw", x=c("rta","iso_o","iso_d"), 
#' vce_robust=TRUE, data=Gravity_zeros)
#' 
#' # Example for data without zero trade flows
#' data(Gravity_no_zeros)
#' 
#' Gravity_no_zeros$lgdp_o <- log(Gravity_no_zeros$gdp_o)
#' Gravity_no_zeros$lgdp_d <- log(Gravity_no_zeros$gdp_d)
#' 
#' NBPML(y="flow", dist="distw", x=c("rta","lgdp_o","lgdp_d"), 
#' vce_robust=TRUE, data=Gravity_no_zeros)
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
#' NBPML(y="flow", dist="distw", x=c("rta","lgdp_o","lgdp_d"), vce_robust=TRUE, data=grav_small_zeros)
#' }
#' 
#' @return
#' The function returns the summary of the estimated gravity model similar to a
#' \code{\link[stats]{glm}}-object.
#' 
#' @seealso \code{\link[MASS]{glm.nb}}, \code{\link[lmtest]{coeftest}}, 
#' \code{\link[sandwich]{vcovHC}}
#' 
#' @export 
#' 
NBPML <- function(y, dist, x, vce_robust=TRUE, data, ...){
  
  if(!is.data.frame(data))                                                stop("'data' must be a 'data.frame'")
  if((vce_robust %in% c(TRUE, FALSE)) == FALSE)                           stop("'vce_robust' has to be either 'TRUE' or 'FALSE'")
  if(!is.character(y)     | !y%in%colnames(data)     | length(y)!=1)      stop("'y' must be a character of length 1 and a colname of 'data'")
  if(!is.character(dist)  | !dist%in%colnames(data)  | length(dist)!=1)   stop("'dist' must be a character of length 1 and a colname of 'data'")
  if(!is.character(x)     | !all(x%in%colnames(data)))                    stop("'x' must be a character vector and all x's have to be colnames of 'data'")  
 
  # Transforming data, logging flows, distances --------------------------------
  
  d                  <- data
  d$dist_log         <- (log(d[dist][,1]))
  d$y                <- d[y][,1] 
  
  # Model ----------------------------------------------------------------------
  
  vars               <- paste(c("dist_log", x), collapse=" + ")
  form               <- paste("y", "~",vars)
  form2              <- stats::as.formula(form)
  model.NBPML        <- MASS::glm.nb(form2, data = d, link = "log")
  model.NBPML.robust <- lmtest::coeftest(model.NBPML, vcov=sandwich::vcovHC(model.NBPML, "HC1"))
  
  # Return --------------------------------------------------------------------- 
  
  if(vce_robust == TRUE){
    summary.NBPML.1                <- .robustsummary.lm(model.NBPML, robust=TRUE)
    summary.NBPML.1$coefficients   <- model.NBPML.robust[1:length(rownames(model.NBPML.robust)),]
    return.object.1                <- summary.NBPML.1
    return.object.1$call           <- form2
    return.object.1$r.squared      <- NULL 
    return.object.1$adj.r.squared  <- NULL
    return.object.1$fstatistic     <- NULL
    return(return.object.1)
  }
  
  if(vce_robust == FALSE){
    return.object.1                <- summary(model.NBPML)
    return.object.1$call           <- form2
    return(return.object.1)}
}