\name{em.cens}
\alias{em.cens}
\title{ Fit univariate censored linear regression model with Normal or Student-t errors }
\description{
Return EM algorithm output for censored regression model (under Normal or Student-t distribution) and calculate diagnostic measures in order to detect influential observations.
}
\usage{em.cens(cc, x, y, nu, type, diagnostic, typediag)}
\arguments{
  \item{cc}{ Vector of censoring indicators. For each observation: 0 if non-censored, 1 if censored }
  \item{x}{ Matrix or vector of covariables }
  \item{y}{ Vector of responde variable }
  \item{nu}{ Initial value for the degree of freedon in case of Student-t model. Must be a positive integer. If type="Normal" nu must not be provided }
  \item{type}{ Distribution to be used. If type="Normal", fit normal model. If type="T", fit Student-t model }
  \item{diagnostic}{ TRUE or FALSE. Indicates if diagnostic measures should or not be computed }
  \item{typediag}{ If diagnostic=TRUE, typediag indicates what kind of diagnostic measure should be computed. If typediag=1, computes generalized Cook distance (and its decomposition into the generalized Cook distance for the parameter subsets: betas and error variance) as a case delection measure. For local influence with case-weight perturbation, set typediag=2. For local influence with scale perturbation, set typediag=3 }
}
\details{ For the Student-t model the specification of the initial value for nu must be made carefully: if the data have many outliers, then you must choose a
small value for nu (near to 2), otherwise you can choose a greater value }
\value{ 
  \item{beta}{ EM estimatives for the beta parameters. The first value correspond to the intercept and the others are relative to each covariable in the order that they appear in the matrix of covariables x }
  \item{sigma2}{ EM estimative for the error variance }
  \item{nu}{ EM estimative for degree of freedom. Only returned when type="T" }
  \item{logver}{ Value of the log-likelihood under the fitted model }
  \item{dp}{ Standard deviation of the estimators of beta }
  \item{measure}{ Gives the vector with calculated values of the diagnostic measure chosen in typediag. Only returned when diagnostic=TRUE }
}
\references{
 Monique B. Massuia, Celso R. Cabral, Larissa A. Matos, Victor H. Lachos. "Influence Diagnostics for Student-t Censored Linear Regression Models".
 }
\author{ Monique Bettio Massuia \email{moniquemassuia@gmail.com}, Larrisa Avila Matos \email{larissaamatos@ime.unicamp.br} and Victor Hugo Lachos \email{hlachos@ime.unicamp.br} }
\examples{ ##see examples in \code{\link{motorettes}} }
\keyword{ censored linear regression }
