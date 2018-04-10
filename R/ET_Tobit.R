#' @title Eaton and Tamura (1994) threshold Tobit model, ET_Tobit
#' 
#' @description \code{ET_Tobit} estimates gravity models in their additive form
#' by conducting a left-censored regression.
#' It follows the Eaton and Tamura (1994) Tobit model,
#' also called threshold Tobit model, where,
#' instead of adding number \code{1} to the dependent variable as done 
#' in \code{\link[gravity]{Tobit}}, the constant added to the
#' data is estimated and interpreted as a threshold.
#' For estimating this threshold, we follow Carson and Sun (2007).
#' 
#' @details \code{ET_Tobit} represents the Eaton and Tamura (1994) Tobit model
#' which is often used when several gravity models are compared.
#' When taking the log of the gravity equation flows equal to zero
#' constitute a problem as their log is not defined.
#' Therefore, a constant is added to the flows. 
#' This constant, opposed to \code{\link[gravity]{Tobit}}, is estimated. 
#' Compared to the usual ET-Tobit approaches, in this package, the estimation
#' of the threshold is done before the other parameters are estimated.
#' We follow Carson and Sun (2007), who show that taking the minimum
#' positive flow value as an estimate of the threshold is super-consistent and that
#' using this threshold estimate ensures that the parameter MLE 
#' are asymptotically normal with the asymptotic variance
#' identical to the variance achieved when the threshold is known. 
#' Hence, first the threshold is estimated as the minimum positive flow. 
#' This threshold is added to the flow variable. It is logged
#' afterwards and taken as the dependent variable. 
#' The Tobit estimation is then conducted using the 
#' \code{\link[censReg]{censReg}} function and setting the lower bound 
#' equal to the log of the minimum positive flow value which was added to all
#' observations.
#' A Tobit regression represents a combination of a binary and a
#' linear regression. 
#' This procedure has to be taken into consideration when
#' interpreting the estimated coefficients.
#' The marginal effects of an explanatory variable on the expected value of 
#' the dependent variable equals the product of both the probability of the 
#' latent variable exceeding the threshold and the marginal effect of the 
#' explanatory variable of the expected value of the latent variable. 
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
#' Up to now, the function is designed for cross-sectional data,
#' but can be easily extended to panel data using the 
#' \code{\link[censReg]{censReg}} function.
#' A robust estimations is not implemented to the present
#' as the \code{\link[censReg]{censReg}} function is not
#' compatible with the \code{\link[sandwich]{vcovHC}} function.
#' 
#' For a more elaborate Tobit function, see \code{\link[gravity]{EK_Tobit}} 
#' for the
#' Eaton and Kortum (2001) Tobit model where each zero trade volume
#' is assigned a country specific interval with the upper 
#' bound equal to the minimum positive trade level of the respective
#' importing country.
#' 
#' @param y name (type: character) of the dependent variable in the dataset 
#' \code{data}, e.g. trade flows.
#' Following Carson and Sun (2007), the smallest positive flow value is 
#' used as an estimate of the threshold. 
#' It is added to \code{y}, the transformed variable is 
#' logged and taken as the dependent variable in the Tobit estimation with 
#' lower bound equal to the log of the smallest possible flow value. 
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
#' @param ... additional arguments to be passed to \code{ET_Tobit}.
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
#' ET_Tobit(y="flow", dist="distw", x=c("rta","lgdp_o","lgdp_d"), 
#' data=Gravity_zeros)
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
#' ET_Tobit(y="flow", dist="distw", x=c("rta","lgdp_o","lgdp_d"), data=grav_small_zeros)
#' }
#' 
#' @return
#' The function returns the summary of the estimated gravity model as a 
#' \code{\link[censReg]{censReg}}-object.
#' 
#' @seealso \code{\link[censReg]{censReg}}
#' 
#' 
#' @export 
#' 
ET_Tobit <- function(y, dist, x, data, ...){
  
  if(!is.data.frame(data))                                                stop("'data' must be a 'data.frame'")
  if(!is.character(y)     | !y%in%colnames(data)     | length(y)!=1)      stop("'y' must be a character of length 1 and a colname of 'data'")
  if(!is.character(dist)  | !dist%in%colnames(data)  | length(dist)!=1)   stop("'dist' must be a character of length 1 and a colname of 'data'")
  if(!is.character(x)     | !all(x%in%colnames(data)))                    stop("'x' must be a character vector and all x's have to be colnames of 'data'")  
  
  # Transforming data, logging flows, distances --------------------------------
  
  # determine minimum positive flow value
  d                          <- data
  d$dist_log                 <- (log(d[dist][,1]))
  flow_min                   <- min(d[y][,1][d[y][,1] != 0])
  flow_min_log               <- log(flow_min)
  d$y_plus_flow_min          <- d[y][,1] + flow_min
  d$y_plus_flow_min_log      <- log(d$y_plus_flow_min)
  
  # Model ----------------------------------------------------------------------
  
  vars                       <- paste(c("dist_log", x), collapse = " + ")
  form                       <- paste("y_plus_flow_min_log", "~", vars)
  form2                      <- stats::as.formula(form)
  model.ET_Tobit             <- censReg::censReg(formula = form2, left = flow_min_log, right = Inf, data = d, start = NULL)
  
  # Return --------------------------------------------------------------------- 
  
  return.object.1            <- summary(model.ET_Tobit)
  return.object.1$call       <- form2
  return(return.object.1)
}

