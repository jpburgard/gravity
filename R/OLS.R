#' @title OLS
#' 
#' @description \code{OLS} estimates gravity models in their traditional form
#' via Ordinary Least Squares (OLS). It does not consider Multilateral 
#' Resistance terms.
#' 
#' @details \code{OLS} estimates gravity models in their traditional, additive, 
#' form via Ordinary Least Squares using the \code{lm} function. 
#' Multilateral Resistance terms are not considered by this function. 
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
#' As the coefficients for the country's incomes were often found to be close to 
#' unitary and unitary income elasticities are in line with some theoretical 
#' foundations on international trade, it is sometimes assumed that the income 
#' elasticities are equal to unity. In order to allow for the estimation with 
#' and without the assumption of unitary income elasticities, the option 
#' \code{uie} is built into \code{OLS} with the default set to  \code{FALSE}. 
#' 
#' \code{OLS} estimation can be used for both, cross-sectional and 
#' panel data. Nonetheless, the function is designed to be consistent with the 
#' Stata code for cross-sectional data provided at the website
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#' The function \code{OLS} was therefore tested for cross-sectional data.
#' For the use with panel data no tests were performed. 
#' Therefore, it is up to the user to ensure that the functions can be applied 
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
#' \code{data}, e.g. trade flows. This variable is logged and taken as the 
#' dependent variable in the estimation.
#' If \code{uie=TRUE} the dependent variable is divided by the product of
#' unilateral incomes \code{inc_o} and \code{inc_d}, e.g. GDPs or GNPs of the 
#' countries of interest and logged afterwards.
#' If \code{uie=FALSE} the dependent variable is logged directly. 
#' The transformed variable is then used as the dependent variable and 
#' the logged income variables are used as independent variables in the 
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
#' Interaction terms can be added.
#' 
#' @param inc_o variable name (type: character) of the income of the country of 
#' origin in the dataset \code{data}. If \code{uie=TRUE}, the dependent variable 
#' \code{y} is divided by the product of the incomes \code{inc_d} and \code{inc_o}.
#' If \code{uie=FALSE}, the incomes are logged and taken as independent 
#' variables in the estimation. 
#' If one wants to use more than one unilateral variable, 
#' e.g. GDP and population, those variables have to be merged into one 
#' variable, e.g. GDP per capita, which can be inserted into \code{inc_o}.
#' 
#' @param inc_d variable name (type: character) of the income of the country of 
#' destination in the dataset \code{data}. If \code{uie=TRUE}, the dependent variable 
#' \code{y} is divided by the product of the incomes \code{inc_d} and \code{inc_o}.
#' If \code{uie=FALSE}, the incomes are logged and taken as independent 
#' variables in the estimation. 
#' If one wants to use more than one unilateral variable, 
#' e.g. GDP and population, those variables have to be merged into one 
#' variable, e.g. GDP per capita, which can be inserted into \code{inc_d}.
#' 
#' @param uie Unitary Income Elasticities (type: logic) determines whether the 
#' parameters are to be estimated assuming unitary income elasticities. 
#' The default value is set to \code{FALSE}. If \code{uie} is set \code{TRUE}, 
#' the flows in the dependent variable \code{y} are divided by the product of 
#' the country pairs' incomes before the estimation. If \code{uie} is set to
#' \code{FALSE}, the income variables are logged and taken as independent
#' variables in the estimation. The variable names for the 
#' incomes should be inserted into \code{inc_o} for the country of origin 
#' and into \code{inc_d} for destination country. 
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
#' When using panel data, a variable for the time may be included in the 
#' dataset. Note that the variable for the time dimension should be of 
#' type: factor. See the references for more information on panel data.
#' 
#' @param ... additional arguments to be passed to \code{OLS}.
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
#' data(Gravity_no_zeros)
#' 
#' OLS(y="flow", dist="distw", x=c("rta", "contig", "comcur"), 
#' inc_o="gdp_o", inc_d="gdp_d", uie=FALSE, 
#' vce_robust=TRUE, data=Gravity_no_zeros)
#' 
#' OLS(y="flow", dist="distw", x=c("rta", "comcur", "contig"), 
#' inc_o="gdp_o", inc_d="gdp_d", uie=TRUE, 
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
#' OLS(y="flow", dist="distw", x=c("rta"), inc_o="gdp_o", inc_d="gdp_d", uie=FALSE, vce_robust=TRUE, data=grav_small)
#' }
#' 
#' @return
#' The function returns the summary of the estimated gravity model as an 
#' \code{\link[stats]{lm}}-object.
#' 
#' @seealso \code{\link[stats]{lm}}, \code{\link[lmtest]{coeftest}}, 
#' \code{\link[sandwich]{vcovHC}}
#' 
#'
#' @export 
#' 

OLS <- function(y, dist, x, inc_d, inc_o, uie=FALSE, vce_robust=TRUE, 
                data, ...){

  if(!is.data.frame(data))                                                stop("'data' must be a 'data.frame'")
  if((vce_robust %in% c(TRUE, FALSE)) == FALSE)                           stop("'vce_robust' has to be either 'TRUE' or 'FALSE'")
  if(!is.character(y)     | !y%in%colnames(data)     | length(y)!=1)      stop("'y' must be a character of length 1 and a colname of 'data'")
  if(!is.character(dist)  | !dist%in%colnames(data)  | length(dist)!=1)   stop("'dist' must be a character of length 1 and a colname of 'data'")
  if(!is.character(x)     | !all(x%in%colnames(data)))                    stop("'x' must be a character vector and all x's have to be colnames of 'data'")  

  if(!is.character(inc_d) | !inc_d%in%colnames(data) | length(inc_d)!=1)  stop("'inc_d' must be a character of length 1 and a colname of 'data'")
  if(!is.character(inc_o) | !inc_o%in%colnames(data) | length(inc_o)!=1)  stop("'inc_o' must be a character of length 1 and a colname of 'data'")
  if((uie %in% c(TRUE, FALSE)) == FALSE)                                  stop("'uie' has to be either 'TRUE' or 'FALSE'")
  
  # Transforming data, logging distances ---------------------------------------
  
  d                  <- data
  d$dist_log         <- (log(d[dist][,1]))
  
  # uie == TRUE ----------------------------------------------------------------
  
  if(uie==TRUE){
    
    # Transforming data --------------------------------------------------------
    
    d$y_inc          <- d[y][,1] / (d[inc_o][,1] * d[inc_d][,1])
    d$y_inc_log      <- log(d$y_inc)
    
    # Model --------------------------------------------------------------------
    
    vars             <- paste(c("dist_log", x), collapse = " + ")
    form             <- paste("y_inc_log", "~", vars)
    form2            <- stats::as.formula(form)
    
    model.OLS        <- stats::lm(form2, data = d)
    model.OLS.robust <- lmtest::coeftest(model.OLS, vcov=sandwich::vcovHC(model.OLS, "HC1"))
  }
  
  # uie == FALSE ---------------------------------------------------------------
  
  if(uie==FALSE){
    
    # Transforming data
    d$y_log         <- log(d[y][,1])
    d$inc_o_log     <- log(d[inc_o][,1])
    d$inc_d_log     <- log(d[inc_d][,1])
    
    # Model --------------------------------------------------------------------
    vars            <- paste(c("dist_log", x), collapse = " + ")
    form            <- paste("y_log", "~", vars, "+ inc_o_log + inc_d_log")
    form2           <- stats::as.formula(form)
    
    model.OLS       <- stats::lm(form2, data = d)
  }
  
  # Else -----------------------------------------------------------------------
  
  if((uie %in% c(TRUE, FALSE))==FALSE){
    stop("uie has to be either TRUE or FALSE")}
  
  # Return ---------------------------------------------------------------------
  
  if(vce_robust == TRUE){
    return.object.1      <- .robustsummary.lm(model.OLS, robust=TRUE)
    return.object.1$call <- form2
    return(return.object.1)}
  
  if(vce_robust == FALSE){
    return.object.1      <- .robustsummary.lm(model.OLS, robust=FALSE)
    return.object.1$call <- form2
    return(return.object.1)}
}