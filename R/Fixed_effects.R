#' @title Fixed_Effects
#' 
#' @description \code{Fixed_Effects} estimates gravity models via
#' OLS and fixed effects for the countries of origin and destination. 
#' These effects catch country specific effects.
#' 
#' @details To account for MR terms, Feenstra (2002) and Feenstra (2004) propose to use 
#' importer and exporter fixed effects. Due to the use of these effects, all 
#' unilateral influences such as GDPs can no longer be estimated. 
#' A disadvantage of the use of \code{Fixed_Effects} is that, when applied to 
#' panel data, the number of country-year or country-pair fixed effects can be 
#' too high for estimation. In addition, no comparative statistics are 
#' possible with \code{Fixed_Effects} as the MR terms are not estimated 
#' explicitly. Nevertheless, Head and Mayer (2014) highlight the importance of 
#' the use of fixed effects. 
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
#' Country specific fixed effects are considered by incorporating
#' \code{"iso_o"} and \code{"iso_d"} in \code{fe}. 
#' By including country specific fixed effects, all monadic effects
#' are captured, including Multilateral Resistance terms. 
#' Therefore, no other unilateral variables such as GDP can be
#' included as independent variables in the estimation.
#' 
#' \code{Fixed_Effects} estimation can be used for both, cross-sectional as well as 
#' panel data. 
#' Nonetheless, the function is designed to be consistent with the 
#' Stata code for cross-sectional data provided at the website
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#' The function \code{Fixed_Effects} was therefore tested for 
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
#' 
#' @param dist name (type: character) of the distance variable in the dataset 
#' \code{data} containing a measure of distance between all pairs of bilateral
#' partners. It is logged automatically when the function is executed. 
#' 
#' @param fe vector of names (type: character) of fixed effects.
#' The default is set to the unilateral identifiers
#' \code{"iso_o"} and \code{"iso_d"} for cross-sectional data. 
#' When using panel data, interaction terms of the iso-codes and time
#' may be added in either \code{fe} or \code{x}. 
#' 
#' @param x vector of names (type: character) of those bilateral variables in 
#' the dataset \code{data} that should be taken as the independent variables 
#' in the estimation. If an independent variable is a dummy variable
#' it should be of type numeric (0/1) in the dataset. If an independent variable is 
#' defined as a ratio, it should be logged. 
#' The fixed effects catch all unilateral effects. Therefore, 
#' no other unilateral variables such as GDP can be
#' included as independent variables in the estimation.
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
#' type: factor. 
#' The time variable can be used as a single dependent variable or interaction 
#' term with other variables such as country identifiers by inserting it into 
#' \code{x} or \code{fe}.
#' See the references for more information on panel data.
#' 
#' @param ... additional arguments to be passed to \code{Fixed_Effects}.
#' 
#' @references 
#' For more information on fixed effects as well as informaton on gravity models, 
#' theoretical foundations and suitable estimation methods in general see
#' 
#' Anderson, J. E. (2010) <DOI:10.3386/w16576>
#' 
#' Head, K. and Mayer, T. (2014) <DOI:10.1016/B978-0-444-54314-1.00003-3>
#' 
#' as well as
#' 
#' Anderson, J. E. (1979) <DOI:10.12691/wjssh-2-2-5>
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
#' Fixed_Effects(y="flow", dist="distw", fe=c("iso_o", "iso_d"), 
#' x=c("rta"), vce_robust=TRUE, data=Gravity_no_zeros)
#' 
#' Fixed_Effects(y="flow", dist="distw", fe=c("iso_o", "iso_d"), 
#' x=c("rta", "comcur", "contig"), vce_robust=TRUE, data=Gravity_no_zeros)
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
#' Fixed_Effects(y="flow", dist="distw", fe=c("iso_o", "iso_d"), x=c("rta"), vce_robust=TRUE, data=grav_small)
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
Fixed_Effects <- function(y, dist, fe=c("iso_o", "iso_d"), x, vce_robust=TRUE, data, ...){
  
  if(!is.data.frame(data))                                                stop("'data' must be a 'data.frame'")
  if((vce_robust %in% c(TRUE, FALSE)) == FALSE)                           stop("'vce_robust' has to be either 'TRUE' or 'FALSE'")
  if(!is.character(y)     | !y%in%colnames(data)     | length(y)!=1)      stop("'y' must be a character of length 1 and a colname of 'data'")
  if(!is.character(dist)  | !dist%in%colnames(data)  | length(dist)!=1)   stop("'dist' must be a character of length 1 and a colname of 'data'")
  if(!is.character(x)     | !all(x%in%colnames(data)))                    stop("'x' must be a character vector and all x's have to be colnames of 'data'")  

  if(!is.character(fe) | !all(unique(unlist(strsplit(fe,c("[:]|[*]"))))%in%colnames(data)) | length(fe)<2)  stop("'fe' must be a character vector of length >=2 and all main variables of the fe's have to be colnames of 'data'")
  
  # Transforming data, logging flows and distances -----------------------------
  
  d                      <- data
  d$dist_log             <- (log(d[dist][,1]))
  d$y_log                <- log(d[y][,1])
  
  # Model ----------------------------------------------------------------------
  
  vars                   <- paste(c("dist_log", x, fe), collapse = " + ")
  vars2                  <- paste(vars)
  form                   <- paste("y_log", "~", vars2)
  form2                  <- stats::as.formula(form)
  model.fe               <- stats::lm(form2, data = d)

  # Return ---------------------------------------------------------------------
  
    return.object.1      <- .robustsummary.lm(model.fe, robust=vce_robust)
    return.object.1$call <- form2
    return(return.object.1)
}