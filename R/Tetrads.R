#' @title Tetrads
#' 
#' @description \code{Tetrads} estimates gravity models 
#' by taking the ratio of the ratio of flows.
#' 
#' @details \code{Tetrads} is an estimation method for gravity models
#' developed by Head et al. (2010) (see the references for
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
#' 
#' The function \code{Tetrads} utilizes the multiplicative form of the
#' gravity equation. After choosing a reference importer \code{k} and 
#' exporter \code{ell} one can eliminate importer and exporter fixed effects 
#' by taking the ratio of ratios. Only those exporters trading with the 
#' reference importer and importers trading with the reference exporter will 
#' remain for the estimation. Therefore, reference countries should
#' preferably be countries which trade with every other country in the dataset. 
#' After restircting the data in this way, \code{Tetrads} estimates the gravity 
#' equation in its additive form by OLS.
#' As, by taking the ratio of ratios, all monadic effects diminish, no
#' unilateral variables such as GDP can be included as independent variables.
#' 
#' \code{Tetrads} estimation can be used for both, cross-sectional as well as 
#' panel data. Nonetheless, the function is designed to be consistent with the 
#' Stata code for cross-sectional data provided on the website
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#' The function \code{Tetrads} was therefore tested for cross-sectional data.
#' If Tetrads is used for panel data, the user may has to drop distance as an
#' independent variable as time-invariant effects drop.
#' For applying \code{Tetrads} to panel data see Head, Mayer and Ries (2010).
#' 
#' @param y name (type: character) of the dependent variable in the dataset 
#' \code{data}, e.g. trade flows. It is logged and
#' taken as the dependent variable in the estimation.
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
#' Unilateral effects drop as the ratio of ratios is taken.
#' 
#' @param k reference importing country, default is set to \code{"USA"}.
#' 
#' @param ell reference exporting country, default is set to \code{"JPN"}.
#' 
#' @param multiway_vcov (type: logic) optional; determines whether a function
#' implementing Cameron, Gelbach, & Miller (2011) multi-way clustering of 
#' variance-covariance matrices in the package \code{multiway_vcov} is used
#' for the estimation. In case \code{multiway_vcov=TRUE}, the 
#' \code{cluster.vcov} function is used. The default value is set to \code{TRUE}.
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
#' @param ... additional arguments to be passed to functions used by 
#' \code{Tetrads}.
#' 
#' @references 
#' For information on \code{Tetrads} see
#' 
#' Cameron, A. C., Gelbach, J. B., and Miller, D. L. (2011) <DOI:10.3386/t0327>
#' 
#' Head, K., Mayer, T., & Ries, J. (2010) <DOI:10.1016/j.jinteco.2010.01.002>
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
#' Baier, S. L. and Bergstrand, J. H. (2009) <DOI:10.1016/j.jinteco.2008.10.004>
#' 
#' Baier, S. L. and Bergstrand, J. H. (2010) in Van Bergeijk, P. A., & Brakman, S. (Eds.) (2010) chapter 4 <DOI:10.1111/j.1467-9396.2011.01000.x>
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
#' Tetrads(y="flow", dist="distw", x=c("rta"), k="USA", ell="JPN", 
#' multiway_vcov=TRUE, data=Gravity_no_zeros)
#' 
#' Tetrads(y="flow", dist="distw", x=c("rta", "comcur", "contig"), 
#' k="USA", ell="JPN", multiway_vcov=FALSE, data=Gravity_no_zeros)
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
#' Tetrads(y="flow", dist="distw", x=c("rta"), k=countries_chosen[1], ell=countries_chosen[2], multiway_vcov=FALSE, data=grav_small)
#' }
#' 
#' @return
#' The function returns the summary of the estimated gravity model as an 
#' \code{\link[stats]{lm}}-object.
#' 
#' @seealso \code{\link[stats]{lm}}, \code{\link[lmtest]{coeftest}}, 
#' \code{\link[multiwayvcov]{cluster.vcov}}
#' 
#' @export 
#' 
Tetrads <- function(y, dist, x, k="USA", ell="JPN", multiway_vcov=TRUE, data, ...){
  
  if(!is.data.frame(data))                                                stop("'data' must be a 'data.frame'")
  if(!is.character(y)     | !y%in%colnames(data)     | length(y)!=1)      stop("'y' must be a character of length 1 and a colname of 'data'")
  if(!is.character(dist)  | !dist%in%colnames(data)  | length(dist)!=1)   stop("'dist' must be a character of length 1 and a colname of 'data'")
  if(!is.character(x)     | !all(x%in%colnames(data)))                    stop("'x' must be a character vector and all x's have to be colnames of 'data'")  
 
  if((multiway_vcov %in% c(TRUE, FALSE)) == FALSE)                        stop("'multiway_vcov' has to be either 'TRUE' or 'FALSE'")
  if(!k%in%data$iso_d)                                                    stop("'k' must be in 'data$iso_d'")
  if(!ell%in%data$iso_o)                                                  stop("'ell' must be in 'data$iso_o'")
  
  # Transforming data, log flows and distances ---------------------------------
  
  d           <- data
  d$dist_log  <- (log(d[dist][,1]))
  d$y_log     <- log(d[y][,1])
  
  # Truncating dataset to only those countries which trade with reference
  # importer and exporter ------------------------------------------------------
  
  d_2a        <- d[d$iso_d == k,] # all iso_o which should stay in dataset
  d_2b        <- d_2a$iso_o
  d_2         <- d[d$iso_o %in% d_2b,] 
  d_3a        <- d_2[d_2$iso_o == ell,] # all iso_d which should stay in dataset
  d_3b        <- d_3a$iso_d
  d_3         <- d_2[d_2$iso_d %in% d_3b,] 
  num.ind.var <- length(x) #independent variables apart from distance
  rm(data); rm(d); rm(d_2a);rm(d_2b); rm(d_2); rm(d_3a); rm(d_3b)
  
  # Taking ratios, ratk --------------------------------------------------------
  
  d_3$lXinratk  <- NA
  d_3$ldistratk <- NA
  
  for(i in unique(d_3$iso_o)){
    
    d_3[d_3$iso_o==i,]$lXinratk <- d_3[d_3$iso_o==i,]$y_log - d_3[d_3$iso_o==i & d_3$iso_d==k,]$y_log
    d_3[d_3$iso_o==i,]$ldistratk <- d_3[d_3$iso_o==i,]$dist_log - d_3[d_3$iso_o==i & d_3$iso_d==k,]$dist_log
  }
  
  # Taking ratios, ratk, for the other independent variables -------------------
  
  ind.var.ratk                  <- list(length=num.ind.var+1)
  ind.var.ratk[[num.ind.var+1]] <- d_3$iso_o
  
  for(j in 1:num.ind.var){
    ind.var.ratk[[j]] <- rep(NA, times=nrow(d_3))
  }
  
  for(j in 1:num.ind.var){
    for(i in unique(d_3$iso_o)){
      ind.var.ratk[[j]][ind.var.ratk[[num.ind.var+1]]==i] <- ((d_3[d_3$iso_o==i,])[x[j]])[,1] - (d_3[d_3$iso_o==i & d_3$iso_d==k,])[x[j]][,1]
    }
  }
  
  for(j in 1:length(x)){
    d_3[x[j]]  <-  ind.var.ratk[[j]]
  }
  
  # Taking the ratio of ratios, rat --------------------------------------------
  
  d_3$lXinrat  <- NA
  d_3$ldistrat <- NA
  
  for(i in unique(d_3$iso_d)){
    
    d_3[d_3$iso_d==i,]$lXinrat  <- d_3[d_3$iso_d==i,]$lXinratk - d_3[d_3$iso_d==i & d_3$iso_o==ell,]$lXinratk
    d_3[d_3$iso_d==i,]$ldistrat <- d_3[d_3$iso_d==i,]$ldistratk - d_3[d_3$iso_d==i & d_3$iso_o==ell,]$ldistratk
  }
  
  d_3$y_log_rat    <- d_3$lXinrat
  d_3$dist_log_rat <- d_3$ldistrat
  
  # Taking the ratio of ratios, rat, for the other independent variables -------
  
  ind.var.rat                  <- list(length=num.ind.var+1)
  ind.var.rat[[num.ind.var+1]] <- d_3$iso_d
  
  for(j in 1:num.ind.var){
    ind.var.rat[[j]] <- rep(NA, times=nrow(d_3))
  }
  
  for(j in 1:num.ind.var){
    for(i in unique(d_3$iso_d)){
      
      ind.var.rat[[j]][ind.var.rat[[num.ind.var+1]]==i] <- ((d_3[d_3$iso_d==i,])[x[j]])[,1] - (d_3[d_3$iso_d==i & d_3$iso_o==ell,])[x[j]][,1]
    }
  }
  
  for(j in 1:length(x)){
    d_3[x[j]]  <-  ind.var.rat[[j]]
  }
  
  # Model ----------------------------------------------------------------------
  
  x_rat <- paste0(x,"_rat")
  
  # new row in dataset for independent _rat variable
  for(j in x){
    l        <- which(x == j)
    rat      <- x_rat[l]
    d_3[rat] <- NA
    d_3[rat] <- d_3[x[l]]
  }
  
  vars  <- paste(c("dist_log_rat", x_rat), collapse=" + ")
  form  <- paste("y_log_rat","~",vars)
  form2 <- stats::as.formula(form)
  
  model.Tetrads        <- stats::lm(form2, data = d_3)
  cluster.formula      <- ~ iso_o + iso_d
  model.Tetrads_vcov   <- multiwayvcov::cluster.vcov(model=model.Tetrads, cluster = cluster.formula)
  model.Tetrads.robust <- lmtest::coeftest(x=model.Tetrads, vcov=model.Tetrads_vcov)
  
  # Return ---------------------------------------------------------------------
  
  if(multiway_vcov == TRUE){
    summary.Ted.1              <- .robustsummary.lm(model.Tetrads, robust=TRUE)
    summary.Ted.1$coefficients <- model.Tetrads.robust[1:length(rownames(model.Tetrads.robust)),]
    return.object.1            <- summary.Ted.1
    return.object.1$call       <- form2
    return(return.object.1)}
  
  if(multiway_vcov == FALSE){
    return.object.1            <- .robustsummary.lm(model.Tetrads, robust=FALSE)
    return.object.1$call       <- form2
    return(return.object.1)}
  
}