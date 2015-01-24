\name{em.cens}
\alias{em.cens}
\title{ Fits Univariate Censored Linear Regression Model With Normal or Student-t Errors }
\description{
Returns EM estimates for right censored regression model (under Normal or Student-t distribution) and calculates some diagnostic measures for detecting influential observations
}
\usage{em.cens(cc, x, y, nu="NULL", dist="Normal", diagnostic="FALSE", typediag="NULL")}
\arguments{
  \item{cc}{ Vector of censoring indicators. For each observation: 0 if non-censored, 1 if censored }
  \item{x}{ Design matrix }
  \item{y}{ Vector with the responde variable }
  \item{nu}{ Initial value for the degree of freedon in case of Student-t model (greater than 2) }
  \item{dist}{ Distribution to be used for the errors. "Normal", for normal or "T" for Student-t }
  \item{diagnostic}{ TRUE or FALSE. Indicates if any diagnostic measure should or not be computed }
  \item{typediag}{ If diagnostic=TRUE, typediag indicates which diagnostic measure should be computed. If typediag=1, computes generalized Cook distance (GD) and its decomposition into the generalized Cook distance for the parameter subsets: betas (GDbeta) and sigma2 (GDsigma2). For local influence with case-weight perturbation, set typediag=2. For local influence with scale perturbation, set typediag=3 }
}
\details{ Despite of this function has been built to deal with right censored response variables, one can easily adapt it to work with left censored responses: set -y and -x to obtain the left censored model fit and any diagnostic measure for it. The specification of the initial value for nu must be made carefully: if the data have many outliers, then you must choose a
small value for nu (greater but near to 2), otherwise you can choose a greater value }
\value{ 
  \item{beta}{ EM estimatives for regression coefficients }
  \item{sigma2}{ EM estimative for the error variance }
  \item{nu}{ EM estimative for degree of freedom. Only returned when type="T" }
  \item{logver}{ Value of the log-likelihood under the fitted model }
  \item{SE}{ Standard error for EM estimators }
  \item{measure}{ Vector with the diagnostic measure chosen in typediag. Only returned when diagnostic=TRUE }
  \item{AIC}{ AIC model selection criteria }
  \item{BIC}{ BIC model selection criteria }
  \item{EDC}{ EDC model selection criteria }
}
\references{
 Monique B. Massuia, Celso R. Cabral, Larissa A. Matos, Victor H. Lachos. "Influence Diagnostics for Student-t Censored Linear Regression Models"
 }
\author{ Monique Bettio Massuia \email{moniquemassuia@gmail.com}, Larrisa Avila Matos \email{larissaamatos@ime.unicamp.br} and Victor Hugo Lachos \email{hlachos@ime.unicamp.br} }
\examples{ ##see examples in \code{\link{wage.rates}} }
\keyword{ censored linear regression }
