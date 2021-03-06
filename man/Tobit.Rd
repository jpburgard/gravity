% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Tobit.R
\name{Tobit}
\alias{Tobit}
\title{Left-censored Tobit model with known threshold}
\usage{
Tobit(y, dist, x, added_constant = 1, data, ...)
}
\arguments{
\item{y}{name (type: character) of the dependent variable in the dataset 
\code{data}, e.g. trade flows.
The number \code{1} is added and the transformed variable is logged and 
taken as the dependent variable in the Tobit estimation with lower bound 
equal to \code{0} as \code{log(1)=0} represents the smallest flows 
in the transformed variable.}

\item{dist}{name (type: character) of the distance variable in the dataset 
\code{data} containing a measure of distance between all pairs of bilateral
partners. It is logged automatically when the function is executed.}

\item{x}{vector of names (type: character) of those bilateral variables in 
the dataset \code{data} that should be taken as the independent variables 
in the estimation. If an independent variable is a dummy variable,
it should be of type numeric (0/1) in the dataset. If an independent variable 
is defined as a ratio, it should be logged. 
Unilateral variables such as country dummies or incomes can be added. 
If unilateral metric variables such as GDPs should be used as independent 
variables, those variables have to be logged first and the 
logged variable can be used in \code{x}.
Interaction terms can be added.}

\item{added_constant}{scalar (type: numeric); represents
the constant to be added to the dependent variable. The default value
is \code{1}.
The minimum of \code{log(y+added_constant)} is taken as the 
left boundary in the Tobit model.
In the often used case of \code{added_constant=1}, the 
dependent variable is left-censored at value \code{0} 
as \code{log(1)=0}.}

\item{data}{name of the dataset to be used (type: character). 
To estimate gravity equations, a square gravity dataset including bilateral 
flows defined by the argument \code{y}, ISO-codes of type character 
(called \code{iso_o} for the country of origin and \code{iso_d} for the 
destination country), a distance measure defined by the argument \code{dist} 
and other potential influences given as a vector in \code{x} are required. 
All dummy variables should be of type numeric (0/1). Missing trade flows as 
well as incomplete rows should be excluded from the dataset. 
Zero trade flows are allowed.}

\item{...}{additional arguments to be passed to \code{Tobit}.}
}
\value{
The function returns the summary of the estimated gravity model as a 
\code{\link[censReg]{censReg}}-object.
}
\description{
\code{Tobit} estimates gravity models in their additive form
by conducting a left-censored regression, which, after adding the
constant \code{1} to the dependent variable, utilizes \code{log(1)=0}
as the censoring value.
}
\details{
\code{Tobit} represents the left-censored Tobit (Tobin, 1958) 
approach utilizing a known censoring threshold
which is often used when several gravity models are compared.
When taking the log of the gravity equation flows equal to zero
constitute a problem as their log is not defined.
Therefore, in the execution of the function the number \code{1}
is added to all flows and the \code{log(flows+1)} is
taken as the dependent variable. 
The Tobit estimation is conducted using the \code{\link[censReg]{censReg}}
function and setting the lower bound equal to \code{0} as
\code{log(1)=0} represents the smallest flows in the transformed
variable. 
A Tobit regression represents a combination of a binary and a
linear regression. 
This procedure has to be taken into consideration when
interpreting the estimated coefficients.
The marginal effects of an explanatory variable on the expected value of 
the dependent variable equals the product of both the probability of the 
latent variable exceeding the threshold and the marginal effect of the 
explanatory variable of the expected value of the latent variable.

To execute the function a square gravity dataset with all pairs of 
countries, ISO-codes for the country of origin and destination, a measure of 
distance between the bilateral partners as well as all 
information that should be considered as dependent an independent 
variables is needed. 
Missing bilateral flows as well as incomplete rows should be 
excluded from the dataset.  
Zero trade flows are allowed. 

Up to now, the function is designed for cross-sectional data,
but can be easily extended to panel data using the 
\code{\link[censReg]{censReg}} function.
A robust estimations is not implemented to the present
as the \code{\link[censReg]{censReg}} function is not
compatible with the \code{\link[sandwich]{vcovHC}} function.

For more elaborate Tobit functions, see \code{\link[gravity]{ET_Tobit}}
for the Eaton and Tamura (1994) threshold Tobit model where instead
of simply adding number \code{1} to the data the threshold is 
estimated or \code{\link[gravity]{EK_Tobit}} for the 
Eaton and Kortum (2001) Tobit model where each zero trade volume
is assigned a country specific interval with the upper 
bound equal to the minimum positive trade level of the respective
importing country.
}
\examples{
\dontrun{
# Example for data with zero trade flows
data(Gravity_zeros)

Gravity_zeros$lgdp_o <- log(Gravity_zeros$gdp_o)
Gravity_zeros$lgdp_d <- log(Gravity_zeros$gdp_d)

Tobit(y="flow", dist="distw", x=c("rta","lgdp_o","lgdp_d"), 
added_constant = 1, data=Gravity_zeros)
}

\dontshow{
# examples for CRAN checks:
# executable in < 5 sec together with the examples above
# not shown to users

data(Gravity_zeros)
Gravity_zeros$lgdp_o <- log(Gravity_zeros$gdp_o)
Gravity_zeros$lgdp_d <- log(Gravity_zeros$gdp_d)

# choose exemplarily 10 biggest countries for check data
countries_chosen_zeros <- names(sort(table(Gravity_zeros$iso_o), decreasing = TRUE)[1:10])
grav_small_zeros <- Gravity_zeros[Gravity_zeros$iso_o \%in\% countries_chosen_zeros,]
Tobit(y="flow", dist="distw", x=c("rta","lgdp_o","lgdp_d"), added_constant = 1, data=grav_small_zeros)
}

}
\references{
For more information on gravity models, theoretical foundations and
estimation methods in general see 

Anderson, J. E. (1979) <DOI:10.12691/wjssh-2-2-5>

Anderson, J. E. (2010) <DOI:10.3386/w16576>

Anderson, J. E. and van Wincoop, E. (2003) <DOI:10.3386/w8079> 

Baier, S. L. and Bergstrand, J. H. (2009) <DOI:10.1016/j.jinteco.2008.10.004>

Baier, S. L. and Bergstrand, J. H. (2010) in Van Bergeijk, P. A., & Brakman, S. (Eds.) (2010) chapter 4 <DOI:10.1111/j.1467-9396.2011.01000.x>

Head, K., Mayer, T., & Ries, J. (2010) <DOI:10.1016/j.jinteco.2010.01.002>

Head, K. and Mayer, T. (2014) <DOI:10.1016/B978-0-444-54314-1.00003-3>

Santos-Silva, J. M. C. and Tenreyro, S. (2006) <DOI:10.1162/rest.88.4.641> 

and the citations therein.


Especially for Tobit models see

Tobin, J. (1958) <DOI:10.2307/1907382>

Eaton, J., & Tamura, A. (1994) <DOI:10.3386/w4758>

Eaton, J., & Kortum, S. (2001) <DOI:10.3386/w8070>.


See \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook} for gravity datasets and Stata code for estimating gravity models.
}
\seealso{
\code{\link[censReg]{censReg}}
}
