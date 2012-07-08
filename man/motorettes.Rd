\name{motorettes}
\docType{data}
\alias{motorettes}
\title{ Accelerated life tests motorettes }
\description{ Accelerated life tests on electrical insulation in motorettes with censoring times 
}
\usage{data(motorettes)}
\format{ 
A data frame with 40 observed times of life tests on electrical insulation in motorettes at four different temperatures. 
}

\encoding{latin1}

\source{ 
Monique B. Massuia, Celso R. Cabral, Larissa A. Matos, Victor H. Lachos. "Influence Diagnostics for Student-t Censored Linear Regression Models". 
}

\examples{
\dontrun{
##Load the data
data(motorettes)

##Set the response y and covariate x
y <- motorettes$y
x <- motorettes$x
cc<- motorettes$cc


##Fits a censored Normal model and calculates Generalized Cook Distance
EM<-em.cens(cc,x,y,type="Normal",diagnostic=TRUE,typediag=1)
#Makes graphs of the diagnostic measures
GD<-EM$measure$GD
GDbeta<-EM$measure$GDbeta
GDsigma2<-EM$measure$GDsigma2
par(mfrow=c(1,3))
plot(GD,main="Generalized Cook Distance")
plot(GDbeta, main="Generalized Cook Distance - Influence in Beta")
plot(GDsigma2,main="Generalized Cook Distance - Influence in Sigma2")

## Note that, under the Normal model, observations 8 and 9 seems 
## to be influential for all the parameter subsets estimates

##Fits a Student-t model, does not calculate diagnostic measures
#Set initial value 
nu <- 4
EM<-em.cens(cc,x,y,nu,type="T",diagnostic=FALSE)


}
}
\keyword{datasets}
